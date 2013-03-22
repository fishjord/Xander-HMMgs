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

import edu.msu.cme.rdp.alignment.hmm.jni.HMMER3;
import edu.msu.cme.rdp.alignment.hmm.jni.HMMER3Hit;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import java.io.File;

/**
 *
 * @author fishjord
 */
public class HMMKmerFilter {

    private static class RefKmer {

        int modelPos;
        int refFileIndex;

        @Override
        public boolean equals(Object obj) {
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final RefKmer other = (RefKmer) obj;
            if (this.modelPos != other.modelPos) {
                return false;
            }
            if (this.refFileIndex != other.refFileIndex) {
                return false;
            }
            return true;
        }

        @Override
        public int hashCode() {
            int hash = 3;
            hash = 37 * hash + this.modelPos;
            hash = 37 * hash + this.refFileIndex;
            return hash;
        }
    }

    public static void main(String[] args) throws Exception {
        if(args.length != 3) {
            System.err.println("HMMKmerFilter <hmm_file> <query_file> <k>");
            System.exit(1);
        }

        HMMER3 hmmer = new HMMER3(args[0]);
        SeqReader reader = new SequenceReader(new File(args[1]));
        int k = Integer.valueOf(args[2]);

        Sequence seq;
        while((seq = reader.readNextSequence()) != null) {
            for(Sequence framedSeq : ProteinUtils.getInstance().allFrames(seq)) {
                String seqString = framedSeq.getSeqString();
                String protSeq = ProteinUtils.getInstance().translateToProtein(seqString, true, 11);

                HMMER3Hit[] hits = hmmer.findHits(protSeq);
                if(hits.length > 0) {
                    String leftKmer = seqString.substring(0, k);
                    String rightKmer = seqString.substring(seqString.length() - k);
                    for(HMMER3Hit hit : hits) {
                        double covered = hit.getSeqEnd() - hit.getSeqStart() / (double)(protSeq.length());
                        if(covered > .9) {
                            System.out.println(
                                    framedSeq.getSeqName().replace("_", "\t") + "\t" +
                                    hit.getModelName() + "\t" +
                                    hit.getBits() + "\t" +
                                    hit.getHmmStart() + "\t" +
                                    hit.getHmmEnd() + "\t" +
                                    hit.getSeqStart() + "\t" +
                                    hit.getSeqEnd() + "\t" +
                                    leftKmer + "\t" +
                                    rightKmer + "\t" +
                                    framedSeq.getSeqString());
                        }
                    }
                }
            }
        }

        reader.close();
    }
}
