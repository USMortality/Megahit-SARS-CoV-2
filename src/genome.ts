import * as fs from 'fs'
import seedrandom from 'seedrandom'
import cliProgress from 'cli-progress'

const bar = new cliProgress.SingleBar({}, cliProgress.Presets.shades_classic)

// Config
const N_LEN = 1000
const ALPHABET = ['A', 'T', 'C', 'G']
const READ_LEN_MIN = 100
const READ_LEN_MAX = 150
const N_READS = 10000
const N_ORGANISMS = 2
const GENOME_ID = 'XX000001'
const GENOME_NAME = 'Random Test Genome'
const SRR_ID = 'SRR00000001'

const ERR_RATE = 0.05 // 1 in 500
const P_GENOME_READS = 0.5

// Use a seeded random, to ensure the genome remains identical.
const rnd = seedrandom(GENOME_ID)

const replaceAt = (str: string, index: number, replacement: string) => {
  return str.substring(0, index) + replacement + str.substring(index + replacement.length)
}

const randomSequence = (len: number, isSticky = false) => {
  let seq = ''
  for (let i = 0; i < len; i++) {
    const rnd_idx = Math.floor((isSticky ? rnd() : Math.random()) * 4)
    seq += ALPHABET[rnd_idx]
  }
  return seq
}

const generateGenome = () => {
  let result = [`>${GENOME_ID} ${GENOME_NAME}`]
  result.push(randomSequence(N_LEN, true))
  return result
}

const hasError = () => Math.random() * (1 / ERR_RATE) < 1
const isGenomeRead = () => Math.random() * (1 / P_GENOME_READS) < 1

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
  if (isGenomeRead()) { // Real Genome Read
    const idx_start = Math.floor(Math.random() * genome.length)
    const idx_end = idx_start + read_len
    const read = genome.substring(idx_start, idx_end)
    return maybeScrambleRead(read)
  } else { // Random Noise Read
    return randomSequence(read_len)
  }
}

const repeat = (x, times: number) => {
  let result = ''
  for (let i = 0; i < times; i++) result += x
  return result
}

const generateReads = (genome: string, depth: number, n_reads: number) => {
  const file = `./out/${SRR_ID}.fastq`
  // Remove file if exists.
  fs.unlink(file, () => { })

  bar.start(depth * n_reads, 0);

  for (let d = 0; d < depth; d++) {
    for (let i = 0; i < n_reads; i++) {
      const result: string[] = []
      let read = ''
      while (read.length <= READ_LEN_MIN) read = generateRead(genome)
      result.push(`@SRR00000000.1 ${i + 1} length=${read.length}`)
      result.push(read)
      result.push(`+SRR00000000.1 ${i + 1} length=${read.length}`)
      result.push(repeat('I', read.length))

      fs.appendFileSync(file, result.join('\n') + '\n', err => {
        if (err) console.error(err)
      })
      bar.update((d + 1) * (i + 1))
    }
  }
  bar.stop()
}

console.log('Generating new genome.')
const genome = generateGenome()
console.log('Saving genome...')
fs.writeFile(`./out/${GENOME_ID}.fa`, genome.join('\n'), err => {
  if (err) console.error(err)
})
console.log('Genome saved.')

console.log('Generating reads for genome...')
generateReads(genome[1], N_ORGANISMS, N_READS)

console.log('Finished!')
