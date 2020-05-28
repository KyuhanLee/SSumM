package ssumm;

import java.text.SimpleDateFormat;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectIterator;
import ssumm.SSumM;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;

import static java.lang.Math.sqrt;

public class ReOne extends SSumM {

    
    public ReOne(String datapath){
        super(datapath);
    }

    
    @Override
    public double dataCostPair(int src, int dst, int n){
        Int2IntOpenHashMap source = superGraph[src];
        Int2IntOpenHashMap dest = superGraph[dst];

        int mergeSize = insideSupernode[src].size() + insideSupernode[dst].size();
        long matrixSize = mergeSize;
        int edgeNum;

        if(n==src){
            edgeNum = (source.get(src) & 0x7FFFFFFF) + 2 * (source.get(dst) & 0x7FFFFFFF) + (dest.get(dst) & 0x7FFFFFFF);
            matrixSize *= (matrixSize - 1);
        }
        else{
            edgeNum = (source.get(n) & 0x7FFFFFFF) + (dest.get(n) & 0x7FFFFFFF);
            matrixSize *= insideSupernode[n].size();
        }
        double weight = (double)edgeNum / matrixSize;
        if(src==n){
            edgeNum/=2;
            matrixSize/=2;
        }
        // rareDrop
        if(weight < (double)1/2){
            return -sparseEdgeCost * edgeNum;
        }
        double dense = 2 * log2(numSuperNodes) + log2(maxWeight) + entropy(weight)*matrixSize;
        double sparse = sparseEdgeCost * edgeNum;

        double mergeCost = (dense<=sparse) ? dense : (-sparse);

        return mergeCost;
    }

    // Further sparsification phase
    public void topNDrop(){
        Int2BooleanOpenHashMap tmp;
        IntOpenHashSet nodeSet = new IntOpenHashSet();
        int edgeNum, src, dst;

        int[][] sedges = new int[numSuperEdges][2];
        double[] cost = new double[numSuperEdges];
        int[] idxs = new int[numSuperEdges];
        int[] medians = new int[numSuperEdges];

        int col = 0;
        int check;
        int size;
        int nSize, blockSize;
        int numE;

        for(int n: snList) {
            for (Int2IntOpenHashMap.Entry e : Int2IntMaps.fastIterable(superGraph[n])) {
                if (e.getIntKey() >= n) {
                    numE = e.getIntValue() / ((e.getIntKey() == n) ? 2 : 1);
                    nSize = insideSupernode[n].size();
                    blockSize = nSize * ((e.getIntKey() == n) ? (nSize - 1) : insideSupernode[e.getIntKey()].size());
                    cost[col] = ((double) numE * 2 / blockSize -1) * (double) numE;
                    sedges[col][0] = n;
                    sedges[col][1] = e.getIntKey();
                    idxs[col] = col;
                    col++;
                }
            }
        }

        double summarySize = numNodes * log2(numSuperNodes) + numSuperEdges * (2 * log2(numSuperNodes) + log2(maxWeight));
        int ss = 0, ee = numSuperEdges-1;

        while(summarySize > targetSummarySize) {
            int edgeLft = numSuperEdges - (int)((targetSummarySize - numNodes * log2(numSuperNodes)) / (2 * log2(numSuperNodes) + log2(maxWeight)));
            int _tmp;
            for (int i = ss; i <= ee; i++) {
                medians[i - ss] = idxs[i];
            }
            int chosen = -1;
            int mediansLeft = ee - ss + 1;
            while (mediansLeft > 1) {
                for (int j = 0; j < mediansLeft; j += 5) {
                    int mlen = (j + 5 < mediansLeft) ? 5 : (mediansLeft - j);
                    if (mlen == 5) {
                        if (cost[medians[j]] > cost[medians[j + 1]]) {
                            _tmp = medians[j];
                            medians[j] = medians[j + 1];
                            medians[j + 1] = _tmp;
                        }
                        if (cost[medians[j + 2]] > cost[medians[j + 3]]) {
                            _tmp = medians[j + 2];
                            medians[j + 2] = medians[j + 3];
                            medians[j + 3] = _tmp;
                        }
                        if (cost[medians[j + 1]] < cost[medians[j + 3]]) {
                            medians[j + 3] = medians[j + 4];
                            if (cost[medians[j + 2]] > cost[medians[j + 3]]) {
                                _tmp = medians[j + 2];
                                medians[j + 2] = medians[j + 3];
                                medians[j + 3] = _tmp;
                            }
                        } else {
                            medians[j + 1] = medians[j + 4];
                            if (cost[medians[j]] > cost[medians[j + 1]]) {
                                _tmp = medians[j];
                                medians[j] = medians[j + 1];
                                medians[j + 1] = _tmp;
                            }
                        }
                        if (cost[medians[j + 1]] < cost[medians[j + 3]]) {
                            chosen = (cost[medians[j + 1]] > cost[medians[j + 2]]) ? medians[j + 1] : medians[j + 2];
                        } else {
                            chosen = (cost[medians[j]] > cost[medians[j + 3]]) ? medians[j] : medians[j + 3];
                        }
                    } else if (mlen == 4) {
                        if (cost[medians[j]] > cost[medians[j + 1]]) {
                            _tmp = medians[j];
                            medians[j] = medians[j + 1];
                            medians[j + 1] = _tmp;
                        }
                        if (cost[medians[j + 2]] > cost[medians[j + 3]]) {
                            _tmp = medians[j + 2];
                            medians[j + 2] = medians[j + 3];
                            medians[j + 3] = _tmp;
                        }
                        if (cost[medians[j]] < cost[medians[j + 2]]) {
                            chosen = (cost[medians[j + 1]] < cost[medians[j + 2]]) ? medians[j + 1] : medians[j + 2];
                        } else {
                            chosen = (cost[medians[j]] > cost[medians[j + 3]]) ? medians[j] : medians[j + 3];
                        }
                    } else if (mlen == 3) {
                        if (cost[medians[j]] > cost[medians[j + 1]]) {
                            _tmp = medians[j];
                            medians[j] = medians[j + 1];
                            medians[j + 1] = _tmp;
                        }
                        if (cost[medians[j]] < cost[medians[j + 2]]) {
                            chosen = (cost[medians[j + 1]] < cost[medians[j + 2]]) ? medians[j + 1] : medians[j + 2];
                        } else chosen = medians[j];
                    } else if (mlen == 2) {
                        chosen = (cost[medians[j]] < cost[medians[j + 1]]) ? medians[j] : medians[j + 1];
                    } else {
                        chosen = medians[j];
                    }
                    medians[j / 5] = chosen;
                }
                mediansLeft = (mediansLeft + 4) / 5;
            }
            chosen = medians[0];
            int _ss = ss, _ee = ee, _mid = ee;
            while (_ss <= _mid) {
                if (cost[idxs[_ss]] < cost[chosen]) {
                    _ss++;
                } else if (cost[idxs[_ss]] > cost[chosen]) {
                    _tmp = idxs[_ss];
                    idxs[_ss] = idxs[_mid];
                    idxs[_mid] = idxs[_ee];
                    idxs[_ee] = _tmp;
                    _ee--;
                    _mid--;
                } else {
                    _tmp = idxs[_ss];
                    idxs[_ss] = idxs[_mid];
                    idxs[_mid] = _tmp;
                    _mid--;
                }
            }
            if (edgeLft <= (_ss - ss)) {
                ee = _ss - 1;
            } else if (edgeLft <= (_ee - ss + 1)) {
                for (int i = ss; i < ss + edgeLft; i++) {
                    src = sedges[idxs[i]][0];
                    dst = sedges[idxs[i]][1];
                    superGraph[src].remove(dst);
                    superGraph[dst].remove(src);
                    numSuperEdges--;
                    if (superGraph[src].size() == 0) {
                        if (isolatedId == -1 || isolatedId == src) {
                            isolatedId = src;
                        } else {
                            insideSupernode[isolatedId].addAll(insideSupernode[src]);
                            superGraph[src] = null;
                            insideSupernode[src] = null;
                            removeSupernode(src);
                        }
                    }
                    if (src != dst && superGraph[dst].size() == 0) {
                        if (isolatedId == -1 || isolatedId == dst) {
                            isolatedId = dst;
                        } else {
                            insideSupernode[isolatedId].addAll(insideSupernode[dst]);
                            superGraph[dst] = null;
                            insideSupernode[dst] = null;
                            removeSupernode(dst);
                        }
                    }
                    summarySize = numNodes * log2(numSuperNodes) + numSuperEdges * (2 * log2(numSuperNodes) + log2(maxWeight));
                    if(summarySize <= targetSummarySize) break;
                }
                break;
            } else {
                for (int i = ss; i <= _ee; i++) {
                    src = sedges[idxs[i]][0];
                    dst = sedges[idxs[i]][1];
                    superGraph[src].remove(dst);
                    superGraph[dst].remove(src);
                    numSuperEdges--;
                    if (superGraph[src].size() == 0) {
                        if (isolatedId == -1 || isolatedId == src) {
                            isolatedId = src;
                        } else {
                            insideSupernode[isolatedId].addAll(insideSupernode[src]);
                            superGraph[src] = null;
                            insideSupernode[src] = null;
                            removeSupernode(src);
                        }
                    }
                    if (src != dst && superGraph[dst].size() == 0) {
                        if (isolatedId == -1 || isolatedId == dst) {
                            isolatedId = dst;
                        } else {
                            insideSupernode[isolatedId].addAll(insideSupernode[dst]);
                            superGraph[dst] = null;
                            insideSupernode[dst] = null;
                            removeSupernode(dst);
                        }
                    }
                    summarySize = numNodes * log2(numSuperNodes) + numSuperEdges * (2 * log2(numSuperNodes) + log2(maxWeight));
                    if(summarySize <= targetSummarySize) break;
                }
                if(summarySize <= targetSummarySize) break;
                ss = _ee + 1;
            }
        }
    }

    public double norm(){
        int[] inv = new int[numNodes];
        double err = 0;
        for(int v: snList){
            for(int _u: insideSupernode[v]){
                inv[_u] = v;
            }
        }
        for(long e: edges){
            if(e < 0) continue;
            int _u = (int)(e >> 32);
            int _v = (int)(e & 0x7FFFFFFFL);
            if(_u == _v) continue;
            int v = inv[_v], u = inv[_u];
            int edgeCount = superGraph[u].getOrDefault(v, 0);
            long sz = insideSupernode[u].size();
            sz *= ((u == v) ? (insideSupernode[u].size() - 1) : insideSupernode[v].size());
            double w = edgeCount / (double)sz;
            err += (1-2*w);
        }

        err *= 2;
        for(int v: snList){
            for(Int2IntMap.Entry u: superGraph[v].int2IntEntrySet()){
                int edgeCount = u.getIntValue();
                err += edgeCount;
            }
        }
        err /= numNodes;
        err /= (numNodes - 1);
        return err;
    }
}
