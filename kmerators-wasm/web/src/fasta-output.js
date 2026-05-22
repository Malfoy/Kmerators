export function kmersFromContigs(contigsFa, k) {
  const records = parseFastaRecords(contigsFa);
  const parts = [];
  for (const record of records) {
    if (record.seq.length < k) continue;
    for (let pos = 0; pos + k <= record.seq.length; pos += 1) {
      parts.push(`>${record.id}:kmer_${pos + 1}-${pos + k}\n`);
      parts.push(record.seq.slice(pos, pos + k));
      parts.push("\n");
    }
  }
  return parts.join("");
}

export async function* kmersFromContigTextChunks(textChunks, k, targetChunkChars = 1024 * 1024) {
  if (!Number.isFinite(k) || k < 1) return;

  const encoder = new TextEncoder();
  let id = null;
  let recordIndex = 0;
  let carry = "";
  let nextStart = 1;
  let out = "";

  for await (const rawLine of fastaLines(textChunks)) {
    const line = rawLine.trim();
    if (!line) continue;

    if (line.startsWith(">")) {
      id = line.slice(1).split(/\s+/)[0] || `contig_${recordIndex + 1}`;
      recordIndex += 1;
      carry = "";
      nextStart = 1;
      continue;
    }

    if (id === null) continue;

    const seq = carry + line.replace(/\s+/g, "").toUpperCase();
    for (let pos = 0; pos + k <= seq.length; pos += 1) {
      const start = nextStart++;
      out += `>${id}:kmer_${start}-${start + k - 1}\n${seq.slice(pos, pos + k)}\n`;
      if (out.length >= targetChunkChars) {
        yield encoder.encode(out);
        out = "";
      }
    }
    carry = k > 1 ? seq.slice(Math.max(0, seq.length - k + 1)) : "";
  }

  if (out) {
    yield encoder.encode(out);
  }
}

async function* fastaLines(textChunks) {
  let pending = "";
  for await (const chunk of textChunks) {
    pending += chunk;
    let start = 0;
    for (let index = 0; index < pending.length; index += 1) {
      if (pending.charCodeAt(index) !== 10) continue;
      let line = pending.slice(start, index);
      if (line.endsWith("\r")) line = line.slice(0, -1);
      yield line;
      start = index + 1;
    }
    pending = pending.slice(start);
  }

  if (pending) {
    yield pending.endsWith("\r") ? pending.slice(0, -1) : pending;
  }
}

function parseFastaRecords(text) {
  const records = [];
  let id = null;
  let seq = [];
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line) continue;
    if (line.startsWith(">")) {
      if (id !== null) records.push({ id, seq: seq.join("") });
      id = line.slice(1).split(/\s+/)[0] || `contig_${records.length + 1}`;
      seq = [];
    } else if (id !== null) {
      seq.push(line.replace(/\s+/g, "").toUpperCase());
    }
  }
  if (id !== null) records.push({ id, seq: seq.join("") });
  return records;
}
