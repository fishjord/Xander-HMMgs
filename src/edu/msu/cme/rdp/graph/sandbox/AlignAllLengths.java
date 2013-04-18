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
public class AlignAllLengths {

    public static void main(String[] args) throws Exception {
        if(args.length != 2) {
            System.err.println("USAGE: AlignAllLengths <hmm> <seqfile>");
            System.exit(1);
        }
        File hmmFile = new File(args[0]);
        File seqFile = new File(args[1]);
        ProfileHMM hmm = HMMER3bParser.readModel(hmmFile);

        HMMScorer scorer = new ViterbiScorer(hmm, 5);
        SeqReader reader = new SequenceReader(seqFile);
        Sequence seq;

        while((seq = reader.readNextSequence()) != null) {
            char[] residues = SeqUtils.getUnalignedSeqString(seq.getSeqString()).toCharArray();
            for(int index = 0;index < residues.length;index++) {
                scorer.consume(residues[index]);
                Object[] tmp = scorer.getBestScore(index + 1);
                int k = (Integer)tmp[0];
                double sc = (Double)tmp[1];
                double cc = sc + Math.log(3.0 / (index + 3)) * 2;
                System.out.println(seq.getSeqName() + "\t" + index + "\t" + residues[index] + "\t" + k + "\t" + sc + "\t" + cc + "\t" + (cc - HMMScorer.getNull1(index)) / HMMScorer.ln2);
            }
        }
        reader.close();
    }

}
