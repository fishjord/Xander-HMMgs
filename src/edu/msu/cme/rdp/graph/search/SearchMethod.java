package edu.msu.cme.rdp.graph.search;

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
import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.alignment.hmm.scoring.ForwardScorer;
import edu.msu.cme.rdp.alignment.hmm.scoring.HMMScorer;
import edu.msu.cme.rdp.graph.filter.CodonWalker;
import edu.msu.cme.rdp.graph.search.SearchResult.SearchDirection;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author fishjord
 */
public abstract class SearchMethod {

    public static class PartialResult {

        public PartialResult() {
        }
        public String maxSeq;
        public String alignment;
        public double maxScore;
    }

    public List<SearchResult> search(SearchTarget target) {
        String framedKmer = target.getKmer();
        int frame = target.getFrame();
        /*
         * This never quite worked right...so I tossed out framing the kmer for now
         * This is just a note to remind me I need to do it at some point (maybe)
        if (frame != 0) {
        framedKmer = kmer.substring(frame, kmer.length() - (3 - frame));
        } else {
        framedKmer = kmer;
        }*/

        long leftTime = System.currentTimeMillis();
        List<PartialResult> leftParts = searchInternal(target.getReverseHmm(), target.getReverseHmm().M() - target.getStartState(), framedKmer, target.isProt(), target.getFilter().new LeftCodonFacade(target.getKmer()), false);
        leftTime = System.currentTimeMillis() - leftTime;

        List<SearchResult> ret = new ArrayList();

        for(PartialResult r : leftParts) {
            String nuclSeq = r.maxSeq + framedKmer;
            String alignment = r.alignment + framedKmer.toUpperCase();
            String protSeq = null;
            String scoringSeq = nuclSeq;

            HMMScorer scorer = new ForwardScorer(target.getForwardHmm(), -1);

            if (target.isProt()) {
                scoringSeq = protSeq = ProteinUtils.getInstance().translateToProtein(nuclSeq, true, 11);
            }

            for (char c : scoringSeq.toCharArray()) {
                scorer.consume(c);
            }

            ret.add(new SearchResult(target, target.getKmer(), nuclSeq, alignment, protSeq, SearchDirection.left, target.getReverseHmm().M() - target.getStartState(), r.maxScore, scorer.getMaxScore(), leftTime));
        }

        long rightTime = System.currentTimeMillis();
        List<PartialResult> rightParts = searchInternal(target.getForwardHmm(), target.getStartState(), framedKmer, target.isProt(), target.getFilter().new RightCodonFacade(target.getKmer()), true);
        rightTime = System.currentTimeMillis() - rightTime;

        for(PartialResult r : rightParts) {
            String nuclSeq = framedKmer + r.maxSeq;
            String alignment = framedKmer.toUpperCase() + r.alignment;
            String protSeq = null;
            String scoringSeq = nuclSeq;

            HMMScorer scorer = new ForwardScorer(target.getForwardHmm(), -1);

            if (target.isProt()) {
                scoringSeq = protSeq = ProteinUtils.getInstance().translateToProtein(nuclSeq, true, 11);
            }

            for (char c : scoringSeq.toCharArray()) {
                scorer.consume(c);
            }

            ret.add(new SearchResult(target, target.getKmer(), nuclSeq, alignment, protSeq, SearchDirection.right, target.getStartState(), r.maxScore, scorer.getMaxScore(), rightTime));
        }

        return ret;
    }

    protected abstract List<PartialResult> searchInternal(ProfileHMM hmm, int startingState, String framedKmer, boolean isProt, CodonWalker walker, boolean forward);
}
