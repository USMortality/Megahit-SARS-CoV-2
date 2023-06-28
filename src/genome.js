"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var fs = require("fs");
var seedrandom = require("seedrandom");
var cli_progress_1 = require("cli-progress");
var bar = new cli_progress_1.SingleBar({}, cli_progress_1.Presets.shades_classic);
// Config
var N_LEN = 10;
var ALPHABET = ['A', 'T', 'C', 'G'];
var READ_LEN_MIN = 2;
var READ_LEN_MAX = 3;
var N_READS = 100;
var N_ORGANISMS = 1;
var GENOME_ID = 'XX000000';
var GENOME_NAME = 'Random Test Genome';
var SRR_ID = 'SRR00000000';
var ERR_RATE = 0.05; // 1 in 500
var P_GENOME_READS = 0.5;
// Use a seeded random, to ensure the genome remains identical.
var rnd = seedrandom(GENOME_ID);
var replaceAt = function (str, index, replacement) {
    return str.substring(0, index) + replacement + str.substring(index + replacement.length);
};
var randomSequence = function (len, isSticky) {
    if (isSticky === void 0) { isSticky = false; }
    var seq = '';
    for (var i = 0; i < len; i++) {
        var rnd_idx = Math.floor((isSticky ? rnd() : Math.random()) * 4);
        seq += ALPHABET[rnd_idx];
    }
    return seq;
};
var generateGenome = function () {
    var result = [">".concat(GENOME_ID, " ").concat(GENOME_NAME)];
    result.push(randomSequence(N_LEN, true));
    return result;
};
var hasError = function () { return Math.random() * (1 / ERR_RATE) < 1; };
var isGenomeRead = function () { return Math.random() * (1 / P_GENOME_READS) < 1; };
var maybeScrambleRead = function (read) {
    var new_read = read;
    for (var i = 0; i < new_read.length; i++) {
        if (!hasError())
            continue;
        var rnd_idx = Math.floor(rnd() * 4);
        new_read = replaceAt(new_read, i, ALPHABET[rnd_idx]);
    }
    return new_read;
};
var generateRead = function (genome) {
    var read_var = READ_LEN_MAX - READ_LEN_MIN;
    var read_len = READ_LEN_MIN + Math.floor(Math.random() * read_var);
    if (isGenomeRead()) { // Real Genome Read
        var idx_start = Math.floor(Math.random() * genome.length);
        var idx_end = idx_start + read_len;
        var read = genome.substring(idx_start, idx_end);
        return maybeScrambleRead(read);
    }
    else { // Random Noise Read
        return randomSequence(read_len);
    }
};
var repeat = function (x, times) {
    var result = '';
    for (var i = 0; i < times; i++)
        result += x;
    return result;
};
var generateReads = function (genome, depth, n_reads) {
    var file = "./out/".concat(SRR_ID, ".fastq");
    // Remove file if exists.
    fs.unlink(file, function () { });
    bar.start(depth * n_reads, 0);
    for (var d = 0; d < depth; d++) {
        for (var i = 0; i < n_reads; i++) {
            var result = [];
            var read = '';
            while (read.length <= READ_LEN_MIN)
                read = generateRead(genome);
            result.push("@SRR00000000.1 ".concat(i + 1, " length=").concat(read.length));
            result.push(read);
            result.push("+SRR00000000.1 ".concat(i + 1, " length=").concat(read.length));
            result.push(repeat('I', read.length));
            fs.appendFileSync(file, result.join('\n') + '\n');
            bar.update((d + 1) * (i + 1));
        }
    }
    bar.stop();
};
console.log('Generating new genome.');
var genome = generateGenome();
console.log('Saving genome...');
fs.writeFile("./out/".concat(GENOME_ID, ".fa"), genome.join('\n'), function (err) {
    if (err)
        console.error(err);
});
console.log('Genome saved.');
console.log('Generating reads for genome...');
generateReads(genome[1], N_ORGANISMS, N_READS);
console.log('Finished!');
