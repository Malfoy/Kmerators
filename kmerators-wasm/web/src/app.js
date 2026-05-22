import { Gzip, strToU8 } from "fflate";
import { compress as zstdCompress, init as initZstd } from "@bokuweb/zstd-wasm";
import { Decompress as ZstdDecompress } from "fzstd";
import { kmersFromContigTextChunks } from "./fasta-output.js";

const OUTPUT_CHUNK_BYTES = 1024 * 1024;

const DEFAULT_REFERENCES = {
  transcriptome: [
    {
      kind: "url",
      name: "Homo_sapiens.GRCh38.cdna.all.fa.gz",
      url: "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
    },
    {
      kind: "url",
      name: "Homo_sapiens.GRCh38.ncrna.fa.gz",
      url: "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz",
    },
  ],
  genome: [
    {
      kind: "url",
      name: "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
      url: "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    },
  ],
};

const form = document.querySelector("#run-form");
const runButton = document.querySelector("#run-button");
const resetButton = document.querySelector("#reset-button");
const transcriptomeFile = document.querySelector("#transcriptome-file");
const genomeFile = document.querySelector("#genome-file");
const defaultTranscriptome = document.querySelector("#default-transcriptome");
const defaultGenome = document.querySelector("#default-genome");
const statusTitle = document.querySelector("#status-title");
const statusDetail = document.querySelector("#status-detail");
const elapsedTime = document.querySelector("#elapsed-time");
const metrics = document.querySelector("#metrics");
const downloads = document.querySelector("#downloads");
const downloadButtons = document.querySelector(".download-buttons");
const reportPreview = document.querySelector("#report-preview");

const progressCards = new Map(
  [...document.querySelectorAll("[data-progress]")].map((card) => [card.dataset.progress, card]),
);

const metricEls = {
  queryKmers: document.querySelector("#metric-query-kmers"),
  retained: document.querySelector("#metric-retained"),
  masked: document.querySelector("#metric-masked"),
  genomeHits: document.querySelector("#metric-genome-hits"),
};

let worker = null;
let result = null;
let timer = null;
let startedAt = 0;
let zstdReady = null;

defaultTranscriptome.addEventListener("change", () => {
  transcriptomeFile.disabled = defaultTranscriptome.checked;
});

defaultGenome.addEventListener("change", () => {
  genomeFile.disabled = defaultGenome.checked;
});

form.addEventListener("submit", (event) => {
  event.preventDefault();
  const data = new FormData(form);
  const query = selectedFile(data.get("query"));
  const kmers = parseKmerLengths(data.get("k"));
  if (!kmers.length) {
    setStatus("Invalid k-mer sizes", "Enter one or more positive integers.", true);
    return;
  }
  const minimizerLength = parseOptionalPositiveInt(data.get("minimizerLength"));
  if (minimizerLength != null && kmers.some((k) => minimizerLength > k)) {
    setStatus("Invalid minimizer", "Minimizer length must be no larger than every k-mer size.", true);
    return;
  }
  const payload = {
    sources: {
      query: sourceFromFile(query),
      transcriptome: referenceSources("transcriptome", data),
      genome: referenceSources("genome", data),
    },
    params: {
      k: kmers[0],
      kmers,
      minimizerLength: minimizerLength || 0,
      maxTranscriptome: Number(data.get("maxTranscriptome")),
      maxGenome: Number(data.get("maxGenome")),
      chunkBytes: Number(data.get("chunkMb")) * 1024 * 1024,
    },
  };

  if (!payload.sources.query) {
    setStatus("Missing file", "Select a query FASTA/FASTQ file.", true);
    return;
  }

  startRun(payload);
});

resetButton.addEventListener("click", () => {
  if (worker) {
    worker.terminate();
    worker = null;
  }
  stopTimer();
  form.reset();
  transcriptomeFile.disabled = false;
  genomeFile.disabled = false;
  result = null;
  runButton.disabled = false;
  document.body.classList.remove("is-error");
  setStatus("Ready", "Select files and run the filter.");
  resetProgress();
  metrics.hidden = true;
  downloads.hidden = true;
  reportPreview.hidden = true;
  reportPreview.textContent = "";
});

downloads.addEventListener("click", async (event) => {
  const button = event.target.closest("[data-download]");
  if (!button || !result) return;
  button.disabled = true;
  try {
    await downloadOutput(button.dataset.download, Number(button.dataset.run || 0));
  } catch (error) {
    setStatus("Download error", error?.message || String(error), true);
  } finally {
    button.disabled = false;
  }
});

function startRun(payload) {
  if (worker) worker.terminate();
  result = null;
  runButton.disabled = true;
  resetProgress();
  metrics.hidden = true;
  downloads.hidden = true;
  reportPreview.hidden = true;
  document.body.classList.remove("is-error");
  setStatus("Starting", "Loading WebAssembly worker.");
  startTimer();

  worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });
  worker.addEventListener("message", handleWorkerMessage);
  worker.addEventListener("error", (event) => {
    finishWithError(event.message || "Worker failed.");
  });
  worker.postMessage({ type: "run", payload });
}

function handleWorkerMessage(event) {
  const message = event.data;
  if (message.type === "status") {
    setStatus(message.title, message.detail || "");
    return;
  }
  if (message.type === "progress") {
    setProgress(message.phase, message.loaded, message.total);
    return;
  }
  if (message.type === "result") {
    result = message.result;
    stopTimer();
    runButton.disabled = false;
    setStatus("Done", "Outputs are ready.");
    renderResult(result);
    worker?.terminate();
    worker = null;
    return;
  }
  if (message.type === "error") {
    finishWithError(message.message);
  }
}

function finishWithError(message) {
  stopTimer();
  runButton.disabled = false;
  setStatus("Error", message || "Run failed.", true);
  worker?.terminate();
  worker = null;
}

function renderResult(data) {
  metricEls.queryKmers.textContent = formatCount(data.query_kmers);
  metricEls.retained.textContent = formatCount(data.retained_kmers);
  metricEls.masked.textContent = formatCount(data.masked_kmers);
  metricEls.genomeHits.textContent = formatMaybeCount(data.genome_hits);
  metrics.hidden = false;
  renderDownloadButtons(data);
  downloads.hidden = false;
  reportPreview.textContent = data.files.report_md;
  reportPreview.hidden = false;
}

function setStatus(title, detail, isError = false) {
  statusTitle.textContent = title;
  statusDetail.textContent = detail;
  document.body.classList.toggle("is-error", isError);
}

function resetProgress() {
  for (const phase of progressCards.keys()) {
    setProgress(phase, 0, 0);
  }
}

function setProgress(phase, loaded, total) {
  const card = progressCards.get(phase);
  if (!card) return;
  const output = card.querySelector("output");
  const bar = card.querySelector(".progress-track div");
  const percent = total > 0 ? Math.min(100, (loaded / total) * 100) : loaded > 0 ? 100 : 0;
  output.textContent = total > 0 ? `${formatBytes(loaded)} / ${formatBytes(total)}` : formatBytes(loaded);
  bar.style.width = `${percent}%`;
}

function startTimer() {
  startedAt = performance.now();
  stopTimer();
  timer = setInterval(() => {
    elapsedTime.textContent = `${((performance.now() - startedAt) / 1000).toFixed(1)}s`;
  }, 100);
}

function stopTimer() {
  if (timer) {
    clearInterval(timer);
    timer = null;
  }
}

function renderDownloadButtons(data) {
  const runs = data.runs?.length ? data.runs : [data];
  const buttons = [{ key: "kmers_fa", label: "kmers.fa" }, { key: "contigs_fa", label: "contigs.fa" }, { key: "masked_fa", label: "masked.fa" }];
  downloadButtons.textContent = "";

  if (runs.length === 1) {
    for (const button of buttons) {
      downloadButtons.appendChild(downloadButton(button.key, button.label, 0));
    }
  } else {
    for (const [idx, run] of runs.entries()) {
      const group = document.createElement("div");
      group.className = "download-group";
      const label = document.createElement("span");
      label.textContent = `k=${run.kmer_length}`;
      group.appendChild(label);
      for (const button of buttons) {
        group.appendChild(downloadButton(button.key, button.label, idx));
      }
      downloadButtons.appendChild(group);
    }
  }

  downloadButtons.appendChild(downloadButton("report_md", "report.md", 0));
}

function downloadButton(key, label, runIndex) {
  const button = document.createElement("button");
  button.type = "button";
  button.dataset.download = key;
  button.dataset.run = String(runIndex);
  button.textContent = label;
  return button;
}

async function downloadOutput(key, runIndex = 0) {
  const fileNames = {
    kmers_fa: "kmers.fa",
    contigs_fa: "contigs.fa",
    masked_fa: "masked.fa",
    report_md: "report.md",
  };
  const runs = result.runs?.length ? result.runs : [result];
  const run = runs[runIndex] || runs[0];
  const name = key === "report_md" ? fileNames[key] : outputFileName(run, fileNames[key], runs.length);
  const format = selectedOutputFormat();
  setStatus("Preparing download", name);

  if (key === "report_md") {
    await downloadGeneratedBytes(name, singleChunk(strToU8(result.files.report_md || "")), format, "text/plain;charset=utf-8");
  } else if (key === "contigs_fa") {
    await downloadContigs(run, name, format);
  } else if (key === "kmers_fa") {
    await downloadGeneratedBytes(
      name,
      kmersFromContigTextChunks(textChunksFromBytes(contigByteChunks(run)), run.kmer_length),
      format,
      "text/plain;charset=utf-8",
    );
  } else {
    await downloadGeneratedBytes(
      name,
      singleChunk(strToU8(run.files[key] || "")),
      format,
      "text/plain;charset=utf-8",
    );
  }

  setStatus("Done", "Outputs are ready.");
}

async function downloadContigs(run, name, format) {
  if (format === "zst") {
    downloadBytes(`${name}.zst`, run.contigs_zst, "application/zstd");
    return;
  }

  await downloadGeneratedBytes(name, contigByteChunks(run), format, "text/plain;charset=utf-8");
}

async function downloadGeneratedBytes(name, byteChunks, format, type) {
  if (format === "gz") {
    await downloadBlobFromChunks(`${name}.gz`, gzipChunks(byteChunks), "application/gzip");
    return;
  }

  if (format === "zst") {
    await ensureZstdReady();
    const bytes = await collectBytes(byteChunks);
    downloadBytes(`${name}.zst`, zstdCompress(bytes, -1), "application/zstd");
    return;
  }

  await downloadBlobFromChunks(name, byteChunks, type);
}

async function* contigByteChunks(run) {
  const compressed = run?.contigs_zst;
  if (!compressed) return;

  let pending = [];
  const decompressor = new ZstdDecompress((chunk) => {
    if (chunk?.byteLength) pending.push(chunk);
  });

  for (let offset = 0; offset < compressed.byteLength; offset += OUTPUT_CHUNK_BYTES) {
    const end = Math.min(offset + OUTPUT_CHUNK_BYTES, compressed.byteLength);
    pending = [];
    decompressor.push(compressed.subarray(offset, end), end === compressed.byteLength);
    for (const chunk of pending) {
      yield chunk;
    }
    if (offset && offset % (16 * OUTPUT_CHUNK_BYTES) === 0) {
      await new Promise((resolve) => setTimeout(resolve, 0));
    }
  }
}

async function* textChunksFromBytes(byteChunks) {
  const textDecoder = new TextDecoder();
  for await (const chunk of byteChunks) {
    const text = textDecoder.decode(chunk, { stream: true });
    if (text) yield text;
  }
  const tail = textDecoder.decode();
  if (tail) yield tail;
}

async function* gzipChunks(byteChunks) {
  let pending = [];
  const gzip = new Gzip({ mtime: 0 }, (chunk) => {
    if (chunk?.byteLength) pending.push(chunk);
  });

  for await (const chunk of byteChunks) {
    pending = [];
    gzip.push(chunk, false);
    for (const gzipChunk of pending) {
      yield gzipChunk;
    }
  }

  pending = [];
  gzip.push(new Uint8Array(0), true);
  for (const gzipChunk of pending) {
    yield gzipChunk;
  }
}

async function* singleChunk(bytes) {
  yield bytes;
}

async function collectBytes(byteChunks) {
  const chunks = [];
  let total = 0;
  for await (const chunk of byteChunks) {
    chunks.push(chunk);
    total += chunk.byteLength;
  }

  const bytes = new Uint8Array(total);
  let offset = 0;
  for (const chunk of chunks) {
    bytes.set(chunk, offset);
    offset += chunk.byteLength;
  }
  return bytes;
}

async function downloadBlobFromChunks(name, byteChunks, type) {
  const parts = [];
  for await (const chunk of byteChunks) {
    if (chunk.byteLength) parts.push(chunk);
  }
  downloadBlob(name, new Blob(parts, { type }));
}

function ensureZstdReady() {
  if (!zstdReady) zstdReady = initZstd();
  return zstdReady;
}

function selectedOutputFormat() {
  return document.querySelector('input[name="outputFormat"]:checked')?.value || "plain";
}

function outputFileName(run, name, runCount) {
  return runCount > 1 ? `k${run.kmer_length}.${name}` : name;
}

function downloadBytes(name, bytes, type) {
  downloadBlob(name, new Blob([bytes], { type }));
}

function downloadBlob(name, blob) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = name;
  link.click();
  URL.revokeObjectURL(url);
}

function selectedFile(value) {
  return value instanceof File && value.name ? value : null;
}

function parseKmerLengths(value) {
  const seen = new Set();
  const kmers = [];
  for (const part of String(value || "").split(/[,\s]+/)) {
    if (!part) continue;
    const k = Number(part);
    if (!Number.isInteger(k) || k <= 0 || seen.has(k)) return [];
    seen.add(k);
    kmers.push(k);
  }
  return kmers;
}

function parseOptionalPositiveInt(value) {
  if (value == null || String(value).trim() === "") return null;
  const number = Number(value);
  return Number.isInteger(number) && number > 0 ? number : null;
}

function sourceFromFile(file) {
  return file ? { kind: "file", name: file.name, file } : null;
}

function referenceSources(kind, data) {
  const useDefault = data.get(kind === "genome" ? "defaultGenome" : "defaultTranscriptome") === "on";
  if (useDefault) return DEFAULT_REFERENCES[kind];
  const file = selectedFile(data.get(kind));
  return file ? [sourceFromFile(file)] : [];
}

function formatBytes(bytes) {
  if (!Number.isFinite(bytes) || bytes <= 0) return "0 B";
  const units = ["B", "KB", "MB", "GB"];
  let value = bytes;
  let unit = 0;
  while (value >= 1024 && unit < units.length - 1) {
    value /= 1024;
    unit += 1;
  }
  return `${value.toFixed(value >= 10 || unit === 0 ? 0 : 1)} ${units[unit]}`;
}

function formatCount(value) {
  return new Intl.NumberFormat().format(value || 0);
}

function formatMaybeCount(value) {
  return value == null ? "NA" : formatCount(value);
}
