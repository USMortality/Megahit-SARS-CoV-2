import * as fs from 'fs/promises'
import { SingleBar, Presets } from 'cli-progress'
import { randomSequence, repeat, reverseComplement, reverseString } from './utils'

const bar = new SingleBar({}, Presets.shades_classic)

// Config
const N_READS = 28000000
const SRA_ID = 'SRR00000001'

const generateReads = async () => {
  const file1 = `./out/${SRA_ID}_1.fastq`
  const file2 = `./out/${SRA_ID}_2.fastq`
  // Remove file if exists.
  try { await fs.unlink(file1) } catch (e) { }
  try { await fs.unlink(file2) } catch (e) { }

  bar.start(N_READS, 0);

  for (let i = 0; i < N_READS; i++) {
    const target = randomSequence(Math.round(50 + Math.random() * 100))
    // Forward Read
    const read1 = target
    const result1: string[] = []
    result1.push(`@${SRA_ID}_0:0:0_1:0:0_${i + 1}/1`)
    result1.push(read1)
    result1.push(`+`)
    result1.push(repeat('2', read1.length))
    await fs.appendFile(file1, result1.join('\n') + '\n')

    // Forward Read
    const read2 = reverseString(target)
    const result2: string[] = []
    result2.push(`@${SRA_ID}_0:0:0_1:0:0_${i + 1}/2`)
    result2.push(reverseComplement(read2))
    result2.push(`+`)
    result2.push(repeat('2', read2.length))
    await fs.appendFile(file2, result2.join('\n') + '\n')

    bar.update(i + 1)
  }
  bar.stop()
}

console.log('Generating random reads...')
generateReads()

console.log(`Finished generating Reads: ${SRA_ID}`)
