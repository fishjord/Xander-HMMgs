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
package edu.msu.cme.rdp.graph.utils;

import edu.msu.cme.rdp.alignment.hmm.HMMER3bParser;
import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.alignment.hmm.scoring.ForwardScorer;
import edu.msu.cme.rdp.graph.search.SearchResult.SearchDirection;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.readers.IndexedSeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author fishjord
 */
public class ContigMerger {

    private static class MergedContig implements Comparable<MergedContig> {

        String protSeq;
        String nuclSeq;
        double score;
        String leftContig;
        String rightContig;
        int length;

        public int compareTo(MergedContig o) {
            return Double.compare(o.score, score);
        }
    }

    private static List<MergedContig> mergeAllContigs(Map<String, Sequence> leftContigs, Map<String, Sequence> rightContigs, String kmer, ProfileHMM hmm) {
        List<MergedContig> ret = new ArrayList();
        int k = kmer.length();

        for (Sequence leftContig : leftContigs.values()) {
            String leftSeq = leftContig.getSeqString();
            leftSeq = leftSeq.substring(0, leftSeq.length() - k);

            for (Sequence rightContig : rightContigs.values()) {
                String seq = leftSeq + rightContig.getSeqString();

                MergedContig mergedContig = new MergedContig();
                if (hmm.getAlphabet() == SequenceType.Protein) {
                    mergedContig.nuclSeq = seq;
                    seq = mergedContig.protSeq = ProteinUtils.getInstance().translateToProtein(seq, true, 11);
                } else if (hmm.getAlphabet() == SequenceType.Nucleotide) {
                    mergedContig.nuclSeq = seq;
                } else {
                    throw new IllegalStateException("Cannot handle hmm alpha " + hmm.getAlphabet());
                }
                mergedContig.score = ForwardScorer.scoreSequence(hmm, seq);
                mergedContig.leftContig = leftContig.getSeqName();
                mergedContig.rightContig = rightContig.getSeqName();
                mergedContig.length = seq.length();

                ret.add(mergedContig);
            }
        }

        return ret;
    }

    private static List<MergedContig> mergeContigs(Map<String, Sequence> leftContigs, Map<String, Sequence> rightContigs, String kmer, ProfileHMM hmm) {
        List<MergedContig> finalized = new ArrayList();
        List<MergedContig> candidates = new ArrayList();
        int k = kmer.length();

        Map<String, Set<MergedContig>> contigToMerges = new HashMap();

        for (Sequence leftContig : leftContigs.values()) {
            String leftSeq = leftContig.getSeqString();
            leftSeq = leftSeq.substring(0, leftSeq.length() - k);

            for (Sequence rightContig : rightContigs.values()) {
                String seq = leftSeq + rightContig.getSeqString();

                MergedContig mergedContig = new MergedContig();
                if (hmm.getAlphabet() == SequenceType.Protein) {
                    mergedContig.nuclSeq = seq;
                    seq = mergedContig.protSeq = ProteinUtils.getInstance().translateToProtein(seq, true, 11);
                } else if (hmm.getAlphabet() == SequenceType.Nucleotide) {
                    mergedContig.nuclSeq = seq;
                } else {
                    throw new IllegalStateException("Cannot handle hmm alpha " + hmm.getAlphabet());
                }
                mergedContig.score = ForwardScorer.scoreSequence(hmm, seq);
                mergedContig.leftContig = leftContig.getSeqName();
                mergedContig.rightContig = rightContig.getSeqName();
                mergedContig.length = seq.length();

                candidates.add(mergedContig);

                for (String seqid : new String[]{mergedContig.leftContig, mergedContig.rightContig}) {
                    if (!contigToMerges.containsKey(seqid)) {
                        contigToMerges.put(seqid, new HashSet());
                    }

                    contigToMerges.get(seqid).add(mergedContig);
                }
            }
        }

        Collections.sort(candidates);

        while (candidates.size() > 0) {
            MergedContig mc = candidates.remove(0);

            candidates.removeAll(contigToMerges.get(mc.leftContig));
            candidates.removeAll(contigToMerges.get(mc.rightContig));

            finalized.add(mc);
        }


        return finalized;
    }

    public static void main(String[] args) throws IOException {
        final BufferedReader hmmgsResultReader;
        final IndexedSeqReader nuclContigReader;
        final double minBits;
        final int minProtLength;
        final Options options = new Options();
        final PrintStream out;
        final ProfileHMM hmm;
        final FastaWriter protSeqOut;
        final FastaWriter nuclSeqOut;
        final boolean prot;
        final boolean all;

        options.addOption("a", "all", false, "Generate all combinations for multiple paths, instead of just the best");
        options.addOption("b", "min-bits", true, "Minimum bits score");
        options.addOption("l", "min-length", true, "Minimum length");
        options.addOption("o", "out", true, "Write output to file instead of stdout");

        try {
            CommandLine line = new PosixParser().parse(options, args);

            if (line.hasOption("min-bits")) {
                minBits = Double.valueOf(line.getOptionValue("min-bits"));
            } else {
                minBits = Double.NEGATIVE_INFINITY;
            }

            if (line.hasOption("min-length")) {
                minProtLength = Integer.valueOf(line.getOptionValue("min-length"));
            } else {
                minProtLength = 0;
            }

            if (line.hasOption("out")) {
                out = new PrintStream(line.getOptionValue("out"));
            } else {
                out = System.err;
            }

            all = line.hasOption("all");

            args = line.getArgs();

            if (args.length != 3) {
                throw new Exception("Unexpected number of arguments");
            }

            hmmgsResultReader = new BufferedReader(new FileReader(new File(args[1])));
            nuclContigReader = new IndexedSeqReader(new File(args[2]));
            hmm = HMMER3bParser.readModel(new File(args[0]));

            prot = (hmm.getAlphabet() == SequenceType.Protein);

            if (prot) {
                protSeqOut = new FastaWriter(new File("prot_merged.fasta"));
            } else {
                protSeqOut = null;
            }
            nuclSeqOut = new FastaWriter(new File("nucl_merged.fasta"));

        } catch (Exception e) {
            new HelpFormatter().printHelp("USAGE: ContigMerger [options] <hmm> <hmmgs_file> <nucl_contig>", options);
            System.err.println("Error: " + e.getMessage());
            System.exit(1);
            throw new RuntimeException("I hate you javac");
        }

        String line;
        SearchDirection lastDir = SearchDirection.left; //So this has an assumption built in
        //It depends on hmmgs always outputting left fragments, then right
        //We can't just use the kmer to figure out if we've switched to another starting point
        //because we allow multiple starting model pos, so two different starting
        //positions can have the same starting kmer

        Map<String, Sequence> leftContigs = new HashMap();
        Map<String, Sequence> rightContigs = new HashMap();

        int contigsMerged = 0;
        int writtenMerges = 0;
        long startTime = System.currentTimeMillis();
        String kmer = null;

        while ((line = hmmgsResultReader.readLine()) != null) {
            if (line.startsWith("#")) {
                continue;
            }

            String[] lexemes = line.trim().split("\t");
            if ((lexemes.length != 8 && lexemes.length != 9) || lexemes[0].equals("-")) {
                System.err.println("Skipping line: " + line);
                continue;
            }

            int index = 0;

            String seqid = lexemes[index++];
            kmer = lexemes[index++];
            int modelStart = Integer.valueOf(lexemes[index++]);

            int nuclLength = Integer.valueOf(lexemes[index++]);

            if (lexemes.length == 9) {
                index++;    //prot length column
            }

            SearchDirection dir = SearchDirection.valueOf(lexemes[index++]);

            if (dir != lastDir) {
                if (dir == SearchDirection.left) {
                    List<MergedContig> mergedContigs;

                    if (all) {
                        mergedContigs = mergeAllContigs(leftContigs, rightContigs, kmer, hmm);
                    } else {
                        mergedContigs = mergeContigs(leftContigs, rightContigs, kmer, hmm);
                    }
                    contigsMerged++;

                    for (MergedContig mc : mergedContigs) {
                        String mergedId = mc.leftContig + "_" + mc.rightContig;
                        out.println(mergedId + "\t" + mc.length + "\t" + mc.score);

                        if (mc.score > minBits && mc.length > minProtLength) {
                            if (prot) {
                                protSeqOut.writeSeq(mergedId, mc.protSeq);
                            }
                            nuclSeqOut.writeSeq(mergedId, mc.nuclSeq);

                            writtenMerges++;
                        }
                    }
                    leftContigs.clear();
                    rightContigs.clear();
                }

                lastDir = dir;
            }

            Sequence seq = nuclContigReader.readSeq(seqid);

            if (dir == SearchDirection.left) {
                leftContigs.put(seqid, seq);
            } else if (dir == SearchDirection.right) {
                rightContigs.put(seqid, seq);
            } else {
                throw new IOException("Cannot handle search direction " + dir);
            }
        }

        if (!leftContigs.isEmpty() || !rightContigs.isEmpty()) {
            for (MergedContig mc : mergeContigs(leftContigs, rightContigs, kmer, hmm)) {
                String mergedId = mc.leftContig + "_" + mc.rightContig;
                out.println(mergedId + "\t" + mc.length + "\t" + mc.score);
                contigsMerged++;

                if (mc.score > minBits && mc.length > minProtLength) {
                    if (prot) {
                        protSeqOut.writeSeq(mergedId, mc.protSeq);
                    }
                    nuclSeqOut.writeSeq(mergedId, mc.nuclSeq);
                    writtenMerges++;
                }
            }
        }

        out.close();
        if (prot) {
            protSeqOut.close();
        }
        nuclSeqOut.close();

        System.err.println("Read in " + contigsMerged + " contigs, wrote out " + writtenMerges + " merged contigs in " + ((double) (System.currentTimeMillis() - startTime) / 1000) + "s");
    }
}
