import assert from "node:assert/strict";
import test from "node:test";
import { kmersFromContigTextChunks, kmersFromContigs } from "../src/fasta-output.js";

test("generates kmers from contigs on demand", () => {
  assert.equal(
    kmersFromContigs(">contig1\nACGTA\n", 3),
    ">contig1:kmer_1-3\nACG\n>contig1:kmer_2-4\nCGT\n>contig1:kmer_3-5\nGTA\n",
  );
});

test("skips contigs shorter than k", () => {
  assert.equal(kmersFromContigs(">short\nAC\n", 3), "");
});

test("streams kmers from chunked contigs", async () => {
  const chunks = kmersFromContigTextChunks([">cont", "ig1\nAC", "GTA\n"], 3, 12);
  let output = "";
  for await (const chunk of chunks) {
    output += new TextDecoder().decode(chunk);
  }
  assert.equal(output, ">contig1:kmer_1-3\nACG\n>contig1:kmer_2-4\nCGT\n>contig1:kmer_3-5\nGTA\n");
});
