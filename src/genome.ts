import * as fs from 'fs/promises'
import seedrandom from 'seedrandom'
import { SingleBar, Presets } from 'cli-progress'
import parseArgs from 'minimist'
import {
  ALPHABET,
  fragmentStr, randomSequence, repeat, reverseComplement, reverseString
} from './utils'

const bar = new SingleBar({}, Presets.shades_classic)
const args = parseArgs(process.argv.slice(2))

if (!args.s && !args.sra) throw new Error('must specify SRA via -s or --sra')
const SRA_ID = args.s
const GENOME_ID = args.g || 'XX' + Math.round(Math.random() * 1000000)

// Config
const N_LEN = 29903
const MIN_LEN = 50
const MAX_LEN = 150
const READ_LEN = MAX_LEN - MIN_LEN // 100
const N_ORGANISMS = 10000
const GENOME_NAME = `Random Test Genome: ${GENOME_ID}`
const ERR_RATE = 0.01

// Use a seeded random, to ensure the genome remains identical.
const rnd = seedrandom(GENOME_ID)

const replaceAt = (str: string, index: number, replacement: string) => {
  return str.substring(0, index) + replacement +
    str.substring(index + replacement.length)
}

const generateGenome = () => {
  let result = [`>${GENOME_ID} ${GENOME_NAME}`]
  result.push(randomSequence(N_LEN, rnd()))
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

const replicateString = (str: string, times: number): string[] => {
  const result: string[] = []
  for (let i = 0; i < times; i++) result.push(str)
  return result
}

const generateReads = async (genome: string) => {
  const file1 = `./out/${SRA_ID}_1.fastq`
  const file2 = `./out/${SRA_ID}_2.fastq`
  // Remove file if exists.
  try { await fs.unlink(file1) } catch (e) { }
  try { await fs.unlink(file2) } catch (e) { }

  const genomes = replicateString(genome, N_ORGANISMS)
  const fragmentedGenomes = fragmentStr(genomes, READ_LEN)

  bar.start(fragmentedGenomes.length, 0);

  for (let i = 0; i < fragmentedGenomes.length; i++) {
    const target = fragmentedGenomes[i]
    if (target.length < MIN_LEN || target.length > MAX_LEN) continue
    // Forward Read
    const read1 = maybeScrambleRead(target)
    const result1: string[] = []
    result1.push(`@${SRA_ID}_0:0:0_1:0:0_${i + 1}/1`)
    result1.push(read1)
    result1.push(`+`)
    result1.push(repeat('I', read1.length))
    await fs.appendFile(file1, result1.join('\n') + '\n')

    // Reverse Read
    const read2 = maybeScrambleRead(reverseString(target))
    const result2: string[] = []
    result2.push(`@${SRA_ID}_0:0:0_1:0:0_${i + 1}/2`)
    result2.push(reverseComplement(read2))
    result2.push(`+`)
    result2.push(repeat('I', read2.length))
    await fs.appendFile(file2, result2.join('\n') + '\n')

    bar.update(i + 1)
  }
  bar.stop()
}

let genome: string[]
if (args.g) {
  console.log(`Using genome ${args.g}`)
  const file = (await fs.readFile(args.g, 'utf8')).split('\n')
  genome = [file[0], file.slice(1).join('')]
} else {
  console.log('Generating new genome.')
  genome = generateGenome()
  console.log('Saving genome...')
  await fs.writeFile(`./out/${GENOME_ID}.fa`, genome.join('\n'))
  console.log('Genome saved.')
}

console.log('Generating reads for genome...')
generateReads(genome[1])

console.log(`Finished generating Reads: ${SRA_ID} for Genome: ${GENOME_ID}`)
