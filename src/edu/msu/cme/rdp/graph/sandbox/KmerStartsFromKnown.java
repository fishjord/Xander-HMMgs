/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.graph.sandbox;

import edu.msu.cme.rdp.kmer.io.KmerStart;
import edu.msu.cme.rdp.kmer.io.KmerStartsWriter;
import edu.msu.cme.rdp.kmer.trie.KmerTrie;
import edu.msu.cme.rdp.kmer.trie.ModelPositionKmerGenerator;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.io.File;
import java.util.Arrays;
import java.util.Date;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.CommandLine;

/**
 *
 * @author fishjord
 */
public class KmerStartsFromKnown {

    private static final Options options = new Options();

    static {
        options.addOption("o", "out", true, "Redirect output to file");
        options.addOption("t", "transl-kmer", false, "Translate nucleotide kmers to protein (output protein start points)");
        options.addOption("T", "transl-table", true, "Translation table to use when translating nucleotide to protein sequences");
    }

    ;

    public static void main(String[] args) throws Exception {
        final KmerStartsWriter out;
        final boolean translQuery;
        final int wordSize;
        final int translTable;

        try {
            CommandLine cmdLine = new PosixParser().parse(options, args);
            args = cmdLine.getArgs();

            if (args.length < 2) {
                throw new Exception("Unexpected number of arguments");
            }

            if (cmdLine.hasOption("out")) {
                out = new KmerStartsWriter(cmdLine.getOptionValue("out"));
            } else {
                out = new KmerStartsWriter(System.out);
            }

            if (cmdLine.hasOption("transl-table")) {
                translTable = Integer.valueOf(cmdLine.getOptionValue("transl-table"));
            } else {
                translTable = 11;
            }

            translQuery = cmdLine.hasOption("transl-kmer");
            wordSize = Integer.valueOf(args[0]);

        } catch (Exception e) {
            new HelpFormatter().printHelp("KmerStartsFromKnown <word_size> [name=]<ref_file> ...", options);
            System.err.println(e.getMessage());
            System.exit(1);
            throw new RuntimeException("Stupid jvm");  //While this will never get thrown it is required to make sure javac doesn't get confused about uninitialized variables
        }

        long startTime = System.currentTimeMillis();

        /*
         * if (args.length == 4) { maxThreads = Integer.valueOf(args[3]); } else
         * {
         */

        //}

        System.err.println("Starting kmer mapping at " + new Date());
        System.err.println("*  References:              " + Arrays.asList(args));
        System.err.println("*  Kmer length:             " + wordSize);


        for (int index = 1; index < args.length; index++) {
            String refName;
            String refFileName = args[index];
            if (refFileName.contains("=")) {
                String[] lexemes = refFileName.split("=");
                refName = lexemes[0];
                refFileName = lexemes[1];
            } else {
                String tmpName = new File(refFileName).getName();
                if (tmpName.contains(".")) {
                    refName = tmpName.substring(0, tmpName.lastIndexOf("."));
                } else {
                    refName = tmpName;
                }
            }

            File refFile = new File(refFileName);

            if (SeqUtils.guessSequenceType(refFile) != SequenceType.Nucleotide) {
                throw new Exception("Reference file " + refFile + " contains " + SeqUtils.guessFileFormat(refFile) + " sequences but expected nucleotide sequences");
            }

            SequenceReader seqReader = new SequenceReader(refFile);
            Sequence seq;

            while ((seq = seqReader.readNextSequence()) != null) {
                if (seq.getSeqName().startsWith("#")) {
                    continue;
                }
                ModelPositionKmerGenerator kmers = new ModelPositionKmerGenerator(seq.getSeqString(), wordSize, SequenceType.Nucleotide);

                for (char[] charmer : kmers) {
                    int pos = kmers.getModelPosition() - 1;
                    if (translQuery) {
                        if (pos % 3 != 0) {
                            continue;
                        } else {
                            pos /= 3;
                        }
                    }

                    String kmer = new String(charmer);

                    out.write(new KmerStart(refName,
                            seq.getSeqName(),
                            seq.getSeqName(),
                            kmer,
                            1,
                            pos,
                            translQuery,
                            (translQuery ? ProteinUtils.getInstance().translateToProtein(kmer, true, translTable) : null)));

                }
            }
            seqReader.close();
        }
        out.close();
    }
}
