export const ALPHABET = ['A', 'C', 'G', 'T']
const REV_ALPHABET = ['T', 'G', 'C', 'A']

export const randomSequence = (len: number, rnd: (() => number) | undefined = undefined) => {
  let seq = ''
  for (let i = 0; i < len; i++) {
    const rnd_idx = Math.floor((rnd ? rnd() : Math.random()) * 4)
    seq += ALPHABET[rnd_idx]
  }
  return seq
}

export const getAverageStrLen = (str: string[]) => {
  let result = 0
  for (const s of str) result += s.length
  return Math.floor(result / str.length)
}

export const splitStr = (str: string): string[] => {
  const randomIndex = Math.floor(Math.random() * str.length)
  return [str.slice(0, randomIndex), str.slice(randomIndex)]
}

export const fragmentStr = (str: string[], target_length: number): string[] => {
  if (getAverageStrLen(str) <= target_length) return str

  const result: string[] = []
  // split each piece in two
  for (const s of str) splitStr(s).forEach((el: string) => {
    if (el.length > 0) result.push(el)
  })

  const len = getAverageStrLen(result)
  if (len > target_length) return fragmentStr(result, target_length)
  else return result
}

export const repeat = (x, times: number) => {
  let result = ''
  for (let i = 0; i < times; i++) result += x
  return result
}

export const reverseString = (str) => {
  let newStr = ""
  for (let i = str.length - 1; i >= 0; i--) newStr += str[i]
  return newStr
}

export const reverseComplement = (str: string): string => {
  const result: string[] = []
  for (const s of str) result.push(REV_ALPHABET[ALPHABET.indexOf(s)])
  return result.join('')
}
