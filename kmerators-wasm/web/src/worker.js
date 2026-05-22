import { Gunzip } from "fflate";
import { compress as zstdCompress, init as initZstd } from "@bokuweb/zstd-wasm";
import { Decompress as ZstdDecompress } from "fzstd";
import * as xzModule from "xz-decompress";

const XzReadableStream = xzModule.XzReadableStream || xzModule.default?.XzReadableStream;
const decoder = new TextDecoder();
const encoder = new TextEncoder();
let wasm = null;
let zstdReady = null;

self.addEventListener("message", async (event) => {
  if (event.data?.type !== "run") return;
  try {
    const result = await run(event.data.payload);
    postMessage({ type: "result", result }, transferBuffers(result));
  } catch (error) {
    postMessage({ type: "error", message: error?.message || String(error) });
  }
});

async function run(payload) {
  wasm = await loadWasm();
  const { exports } = wasm.instance;
  const kmers = payload.params.kmers?.length ? payload.params.kmers : [payload.params.k];
  const sessions = kmers.map((k) =>
    createSession(k, payload.params.minimizerLength || 0, payload.params.maxTranscriptome, payload.params.maxGenome),
  );

  try {
    postStatus("Query", payload.sources.query.name);
    const queryBytes = await readWholeFile(payload.sources.query, payload.params.chunkBytes, "query");
    for (const session of sessions) {
      callBytes(exports.kmerators_set_query_fastx || exports.kmerators_set_query_fasta, session, queryBytes);
    }

    const transcriptomeSources = payload.sources.transcriptome || [];
    const genomeSources = payload.sources.genome || [];
    for (const session of sessions) {
      check(exports.kmerators_set_transcriptome_enabled(session, transcriptomeSources.length ? 1 : 0));
      check(exports.kmerators_set_genome_enabled(session, genomeSources.length ? 1 : 0));
    }

    if (transcriptomeSources.length) {
      await scanSources(
        transcriptomeSources,
        payload.params.chunkBytes,
        "transcriptome",
        (bytes) => {
          for (const session of sessions) {
            callBytes(exports.kmerators_scan_transcriptome_chunk, session, bytes);
          }
        },
        () => {
          for (const session of sessions) {
            check(exports.kmerators_finish_transcriptome_source(session));
          }
        },
      );
    } else {
      postProgress("transcriptome", 0, 0);
    }

    if (genomeSources.length) {
      await scanSources(
        genomeSources,
        payload.params.chunkBytes,
        "genome",
        (bytes) => {
          for (const session of sessions) {
            callBytes(exports.kmerators_scan_genome_chunk, session, bytes);
          }
        },
        () => {
          for (const session of sessions) {
            check(exports.kmerators_finish_genome_source(session));
          }
        },
      );
    } else {
      postProgress("genome", 0, 0);
    }

    postStatus("Finalizing", "Building output files.");
    const runs = [];
    for (const session of sessions) {
      check(exports.kmerators_finish(session));
      runs.push(JSON.parse(readResult(session)));
    }
    return compressRunOutputs(combineRuns(runs, payload.params));
  } finally {
    for (const session of sessions) {
      exports.kmerators_session_free(session);
    }
  }
}

function createSession(k, minimizerLength, maxTranscriptome, maxGenome) {
  const { exports } = wasm.instance;
  const session = exports.kmerators_session_new_with_minimizer
    ? exports.kmerators_session_new_with_minimizer(k, minimizerLength, maxTranscriptome, maxGenome)
    : exports.kmerators_session_new(k, maxTranscriptome, maxGenome);
  if (!session) throw new Error(readLastError());
  return session;
}

async function compressRunOutputs(result) {
  await ensureZstdReady();
  for (const run of result.runs) {
    const contigs = run.files?.contigs_fa || "";
    run.contigs_zst = zstdCompress(encoder.encode(contigs), -1);
    delete run.files.contigs_fa;
  }
  return result;
}

function combineRuns(runs, params) {
  const combinedReport = [
    "# kmerators-wasm report",
    "",
    `- k-mer lengths: ${runs.map((run) => run.kmer_length).join(", ")}`,
    "",
    ...runs.flatMap((run) => [
      `## k=${run.kmer_length}`,
      "",
      run.files.report_md.replace(/^# kmerators-wasm report\n\n/, "").trim(),
      "",
    ]),
  ].join("\n");

  return {
    kmer_lengths: runs.map((run) => run.kmer_length),
    kmer_length: runs[0]?.kmer_length || 0,
    query_sequences: runs[0]?.query_sequences || 0,
    query_kmers: sum(runs, "query_kmers"),
    retained_kmers: sum(runs, "retained_kmers"),
    masked_kmers: sum(runs, "masked_kmers"),
    transcriptome_enabled: runs.some((run) => run.transcriptome_enabled),
    genome_enabled: runs.some((run) => run.genome_enabled),
    transcriptome_hits: maybeSum(runs, "transcriptome_hits"),
    genome_hits: maybeSum(runs, "genome_hits"),
    runs,
    files: { report_md: combinedReport },
  };
}

function sum(runs, key) {
  return runs.reduce((total, run) => total + (run[key] || 0), 0);
}

function maybeSum(runs, key) {
  return runs.some((run) => run[key] != null) ? sum(runs, key) : null;
}

function transferBuffers(result) {
  return result.runs.map((run) => run.contigs_zst?.buffer).filter(Boolean);
}

async function ensureZstdReady() {
  if (!zstdReady) zstdReady = initZstd();
  await zstdReady;
}

async function loadWasm() {
  if (wasm) return wasm;
  const response = await fetch(new URL("../wasm/kmerators_wasm_core.wasm", import.meta.url));
  if (!response.ok) {
    throw new Error(`failed to fetch wasm: ${response.status}`);
  }
  const bytes = await response.arrayBuffer();
  const instance = await WebAssembly.instantiate(bytes, {});
  wasm = instance;
  return wasm;
}

async function readWholeFile(source, chunkBytes, phase) {
  const chunks = [];
  let total = 0;
  let lastLoaded = 0;
  let lastTotal = 0;
  for await (const chunk of decodedSourceChunks(source, chunkBytes)) {
    chunks.push(chunk.bytes);
    total += chunk.bytes.byteLength;
    lastLoaded = chunk.loaded;
    lastTotal = chunk.total;
    postProgress(phase, chunk.loaded, chunk.total);
  }

  const out = new Uint8Array(total);
  let offset = 0;
  for (const chunk of chunks) {
    out.set(chunk, offset);
    offset += chunk.byteLength;
  }
  postProgress(phase, lastTotal > 0 ? lastTotal : total || lastLoaded, lastTotal);
  return out;
}

async function scanSources(sources, chunkBytes, phase, scanChunk, finishSource) {
  let baseLoaded = 0;
  const knownTotal = sources.every((source) => source.size)
    ? sources.reduce((sum, source) => sum + source.size, 0)
    : 0;
  for (const source of sources) {
    postStatus(phaseTitle(phase), source.name);
    let sourceLoaded = 0;
    let sourceTotal = source.size || 0;
    for await (const chunk of decodedSourceChunks(source, chunkBytes)) {
      scanChunk(chunk.bytes);
      sourceLoaded = chunk.loaded;
      sourceTotal = chunk.total || sourceTotal;
      const aggregateTotal = knownTotal || (sourceTotal ? baseLoaded + sourceTotal : 0);
      postProgress(phase, baseLoaded + sourceLoaded, aggregateTotal);
    }
    finishSource();
    baseLoaded += sourceTotal || sourceLoaded;
  }
}

async function* decodedSourceChunks(source, chunkBytes) {
  const name = source.name.toLowerCase();
  if (name.endsWith(".gz")) {
    yield* inflateWithPushDecoder(source, chunkBytes, (onData) => new Gunzip(onData));
    return;
  }
  if (name.endsWith(".zst") || name.endsWith(".zstd")) {
    yield* inflateWithPushDecoder(source, chunkBytes, (onData) => new ZstdDecompress(onData));
    return;
  }
  if (name.endsWith(".xz")) {
    yield* xzDecodedChunks(source);
    return;
  }

  yield* rawSourceChunks(source, chunkBytes);
}

async function* inflateWithPushDecoder(source, chunkBytes, createDecoder) {
  let pending = [];
  const decompressor = createDecoder((data) => {
    if (data?.byteLength) pending.push(data);
  });

  for await (const chunk of rawSourceChunks(source, chunkBytes)) {
    pending = [];
    decompressor.push(chunk.bytes, chunk.loaded === chunk.total);
    for (const bytes of pending) {
      yield { bytes, loaded: chunk.loaded, total: chunk.total };
    }
  }
}

async function* xzDecodedChunks(source) {
  if (!XzReadableStream) {
    throw new Error("xz decoder failed to load in this browser bundle.");
  }
  const stream = new XzReadableStream(sourceStream(source));
  const reader = stream.getReader();
  let loaded = 0;
  while (true) {
    const { value, done } = await reader.read();
    if (done) break;
    loaded += value.byteLength;
    yield { bytes: value, loaded, total: 0 };
  }
}

async function* rawSourceChunks(source, chunkBytes) {
  if (source.kind === "url") {
    yield* rawUrlChunks(source);
    return;
  }

  const file = source.file;
  let offset = 0;
  while (offset < file.size) {
    const end = Math.min(file.size, offset + chunkBytes);
    const bytes = new Uint8Array(await file.slice(offset, end).arrayBuffer());
    offset = end;
    yield { bytes, loaded: offset, total: file.size };
  }
}

async function* rawUrlChunks(source) {
  const response = await fetch(source.url);
  if (!response.ok || !response.body) {
    throw new Error(`${source.name}: HTTP ${response.status}`);
  }
  const total = Number(response.headers.get("Content-Length")) || 0;
  source.size = total;
  const reader = response.body.getReader();
  let loaded = 0;
  while (true) {
    const { value, done } = await reader.read();
    if (done) break;
    loaded += value.byteLength;
    yield { bytes: value, loaded, total };
  }
}

function sourceStream(source) {
  if (source.kind === "file") return source.file.stream();
  return new ReadableStream({
    async start(controller) {
      try {
        for await (const chunk of rawUrlChunks(source)) {
          controller.enqueue(chunk.bytes);
        }
        controller.close();
      } catch (error) {
        controller.error(error);
      }
    },
  });
}

function phaseTitle(phase) {
  return phase === "genome" ? "Genome" : phase === "transcriptome" ? "Transcriptome" : "Query";
}

function callBytes(fn, session, bytes) {
  const { exports } = wasm.instance;
  if (!bytes.byteLength) {
    check(fn(session, 0, 0));
    return;
  }
  if (bytes.byteLength > 0xffffffff) {
    throw new Error("single wasm transfer is larger than the wasm32 address space");
  }
  const ptr = exports.kmerators_alloc(bytes.byteLength);
  const offset = wasmPtr(ptr);
  if (!offset) throw new Error("wasm allocation failed");
  assertWasmRange(offset, bytes.byteLength);
  try {
    new Uint8Array(exports.memory.buffer, offset, bytes.byteLength).set(bytes);
    check(fn(session, ptr, bytes.byteLength));
  } finally {
    exports.kmerators_dealloc(ptr, bytes.byteLength);
  }
}

function check(code) {
  if (code !== 0) {
    throw new Error(readLastError());
  }
}

function readResult(session) {
  const { exports } = wasm.instance;
  const ptr = exports.kmerators_result_ptr(session);
  const offset = wasmPtr(ptr);
  const len = wasmU32(exports.kmerators_result_len(session));
  if (!offset || !len) return "";
  assertWasmRange(offset, len);
  const bytes = new Uint8Array(exports.memory.buffer, offset, len).slice();
  return decoder.decode(bytes);
}

function readLastError() {
  const { exports } = wasm.instance;
  const ptr = exports.kmerators_last_error_ptr();
  const offset = wasmPtr(ptr);
  const len = wasmU32(exports.kmerators_last_error_len());
  if (!offset || !len) return "unknown wasm error";
  assertWasmRange(offset, len);
  const bytes = new Uint8Array(exports.memory.buffer, offset, len).slice();
  return decoder.decode(bytes);
}

function wasmPtr(ptr) {
  return ptr >>> 0;
}

function wasmU32(value) {
  return value >>> 0;
}

function assertWasmRange(offset, len) {
  const end = offset + len;
  const size = wasm.instance.exports.memory.buffer.byteLength;
  if (!Number.isSafeInteger(end) || offset > size || end > size) {
    throw new Error(`wasm memory range ${offset}..${end} is outside ${size} bytes`);
  }
}

function postStatus(title, detail) {
  postMessage({ type: "status", title, detail });
}

function postProgress(phase, loaded, total) {
  postMessage({ type: "progress", phase, loaded, total });
}
