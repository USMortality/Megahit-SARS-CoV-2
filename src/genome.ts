import * as fs from 'fs'
import seedrandom from 'seedrandom'

// Config
const N_LEN = 200
const ALPHABET = ['A', 'T', 'C', 'G']
const READ_LEN_MIN = 120
const READ_LEN_MAX = 150
const N_READS = 100
const N_DEPTH = 10
const GENOME_ID = 'XX000001'
const GENOME_NAME = 'Random Test Genome'
const SRR_ID = 'SRR00000001'

const ERR_RATE = 0.01

// Use a seeded random, to ensure the genome remains identical.
const rnd = seedrandom(GENOME_ID)

const replaceAt = (str: string, index: number, replacement: string) => {
  return str.substring(0, index) + replacement + str.substring(index + replacement.length);
}

const generateGenome = () => {
  let result = [`>${GENOME_ID} ${GENOME_NAME}`]
  let genome = ''
  for (let i = 0; i < N_LEN; i++) {
    const rnd_idx = Math.floor(rnd() * 4)
    genome += ALPHABET[rnd_idx]
  }
  result.push(genome)
  return result
}

const hasError = () => Math.random() * (1 / ERR_RATE) < 1

const maybeScrambleRead = (read: string) => {
  let new_read = read
  for (let i = 0; i < new_read.length; i++) {
    if (!hasError()) continue
    const rnd_idx = Math.floor(rnd() * 4)
    new_read = replaceAt(new_read, i, ALPHABET[rnd_idx])
  }
  return new_read
}

const generateRead = (genome: string) => {
  const read_var = READ_LEN_MAX - READ_LEN_MIN
  const read_len = READ_LEN_MIN + Math.floor(Math.random() * read_var)
  const idx_start = Math.floor(Math.random() * genome.length)
  const idx_end = idx_start + read_len
  const read = genome.substring(idx_start, idx_end)
  return maybeScrambleRead(read)
}

const repeat = (x, times: number) => {
  let result = ''
  for (let i = 0; i < times; i++) result += x
  return result
}

const generateReads = (genome: string, n_reads: number, depth: number) => {
  const reads: string[] = []
  for (let j = 0; j < depth; j++) {
    for (let i = 0; i < n_reads; i++) {
      let read = ''
      while (read.length <= READ_LEN_MIN) read = generateRead(genome)
      reads.push(`@SRR00000000.1 ${i + 1} length=${read.length}`)
      reads.push(read)
      reads.push(`+SRR00000000.1 ${i + 1} length=${read.length}`)
      reads.push(repeat('I', read.length))
    }
  }
  return reads
}

console.log('Generating new genome.')
const genome = generateGenome()
console.log('Generating reads for genome.')
const reads = generateReads(genome[1], N_READS, N_DEPTH)

console.log('Saving...')
fs.writeFile(`./out/${GENOME_ID}.fa`, genome.join('\n'), err => {
  if (err) console.error(err);
})
fs.writeFile(`./out/${SRR_ID}.fastq`, reads.join('\n'), err => {
  if (err) console.error(err);
})
console.log('Done!')
