/*
 * Copyright (C) 2013 Jordan Fish <fishjord at msu.edu>
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

import edu.msu.cme.rdp.alignment.hmm.HMMER3bParser;
import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.alignment.hmm.TSC;
import edu.msu.cme.rdp.alignment.hmm.scoring.HMMScorer;
import edu.msu.cme.rdp.alignment.hmm.scoring.ViterbiScorer;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.io.File;

/**
 *
 * @author Jordan Fish <fishjord at msu.edu>
 */
public class AlignAllLengthsManualViterbi {

    public static void main(String[] args) throws Exception {
        if(args.length != 2) {
            System.err.println("USAGE: AlignAllLengths <hmm> <seqfile>");
            System.exit(1);
        }
        File hmmFile = new File(args[0]);
        File seqFile = new File(args[1]);
        ProfileHMM hmm = HMMER3bParser.readModel(hmmFile);

        SeqReader reader = new SequenceReader(seqFile);
        Sequence seq;

        int startingState = 1;
        double ret = 0;

        while((seq = reader.readNextSequence()) != null) {
            char[] residues = SeqUtils.getUnalignedSeqString(seq.getSeqString()).toCharArray();
            for (int index = 0; index < residues.length; index++) {
                double s = hmm.msc(startingState + index, residues[index]) + hmm.tsc(startingState + index - 1, TSC.MM);
                ret += s;


                int k = startingState + index;
                double sc = ret;
                double cc = sc + Math.log(2.0 / (k + 2)) * 2;
                System.out.println(seq.getSeqName() + "\t" + index + "\t" + residues[index] + "\t" +hmm.msc(startingState + index, residues[index]) + "\t" + hmm.tsc(startingState + index, TSC.MM) + "\t" + s + "\t" + k + "\t" + sc + "\t" + cc + "\t" + (cc - HMMScorer.getNull1(k)) / HMMScorer.ln2);
            }
        }
        reader.close();
    }

}
