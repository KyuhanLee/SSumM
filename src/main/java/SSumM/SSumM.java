package ssumm;

import java.text.SimpleDateFormat;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectIterator;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;

import static java.lang.Math.sqrt;

public class SSumM {
    private String dataPath;
    protected double targetSummarySize;

    protected int numNodes;
    protected int numEdges;

    protected static final double log2 = Math.log(2.0);
    protected static final int divThreshold = 500;

    protected Int2IntOpenHashMap node2Idx = new Int2IntOpenHashMap();
    protected IntArrayList deg = new IntArrayList();

    protected long[] edges;
    protected Int2IntOpenHashMap[] superGraph;

    protected IntArrayList snList;
    protected int[] rep;
    protected IntArrayList[] insideSupernode;

    protected double[] dataCost;
    protected int[] modelCost;

    protected long[][] dup = new long[divThreshold][divThreshold];
    protected long cntFlag;
    protected int numSuperNodes;
    protected int numSuperEdges;
    protected int isolatedId = -1;
    protected double sparseEdgeCost;
    protected int maxWeight;

    ThreadLocalRandom rnd = ThreadLocalRandom.current();

    protected SSumM(String dataPath){
        this.dataPath = dataPath;
        this.cntFlag = 0;
    }


    public void addNode(int v){
        int idx = node2Idx.getOrDefault(v, -1);
        if(idx < 0) {
            node2Idx.put(v, numNodes);
            idx = numNodes;
            deg.add(0);
            numNodes++;
        }
        deg.set(idx, deg.getInt(idx) + 1);
    }

    public boolean checkFormat(String number) {
        boolean isInt = false;
        try {
            Integer.parseInt(number);
            isInt = true;
        } catch (NumberFormatException e) {
        }
        return isInt;
    }

    public void inputGraph(){
        String line;
        String[] parts;
        int src, dst;
        int cnt = 0;

        SimpleDateFormat format = new SimpleDateFormat( "yyyy-MM-dd HH:mm:ss");
        String formatTime = format.format(System.currentTimeMillis());
        System.out.println("Data Read Start: " + formatTime);
        try {
            BufferedReader br = new BufferedReader(new FileReader(dataPath));
            while ((line = br.readLine()) != null) {
                parts = line.split("\\s");
                if ((parts.length >= 2) && checkFormat(parts[0]) && checkFormat(parts[1])) {
                    try {
                        src = Integer.parseInt(parts[0]); addNode(src);
                        dst = Integer.parseInt(parts[1]); addNode(dst);
                        cnt++;
                    } catch (NumberFormatException e) {
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            System.err.println(e);
        }
        edges = new long[cnt];
        superGraph = new Int2IntOpenHashMap[numNodes];
        for(int i=0;i<numNodes;i++) superGraph[i] = new Int2IntOpenHashMap(deg.getInt(i));
        insideSupernode = new IntArrayList[numNodes];
        for(int i=0;i<numNodes;i++) insideSupernode[i] = new IntArrayList(new int[]{i});
        deg = null;

        rep = new int[numNodes];
        snList = new IntArrayList(numNodes);
        for(int i=0;i<numNodes;i++){
            rep[i] = i;
            snList.add(i);
        }

        try {
            BufferedReader br = new BufferedReader(new FileReader(dataPath));
            while ((line = br.readLine()) != null) {
                parts = line.split("\\s");
                if ((parts.length >= 2) && checkFormat(parts[0]) && checkFormat(parts[1])) {
                    try {
                        src = node2Idx.get(Integer.parseInt(parts[0]));
                        dst = node2Idx.get(Integer.parseInt(parts[1]));
                        if(src == dst) continue;
                        if(superGraph[src].put(dst, 1) == 0) {
                            edges[numEdges++] = (((long)src) << 32) + dst;
                            superGraph[dst].put(src, 1);
                        }
                    } catch (Exception e) {
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            System.err.println(e);
        }

        System.out.println("|V|: " + numNodes);
        System.out.println("|E|: " + numEdges);
        for(int i=numEdges;i<cnt;i++) edges[i] = -1;

        numSuperNodes = numNodes;
        numSuperEdges = numEdges;

        sparseEdgeCost = 2*log2((double)numNodes);
        maxWeight = 1;
        dataCost = new double[numNodes];
        for(int i=0;i<numNodes;i++) dataCost[i] = 0;
        modelCost = new int[numNodes];
        for(int i=0;i<numNodes;i++) modelCost[i] = superGraph[i].size();
        System.out.println("Finished reading the input graph");
    }

    // Returns minimimum data cost of an edge between src & dst
    public double dataCostPair(int src, int dst){
        double tmpDataCost;

        int tmpEdgeNum = (superGraph[src].get(dst) & 0x7FFFFFFF);

        if(superGraph[src].get(dst) > 0){
            int srcSize = insideSupernode[src].size();
            long matrixSize = srcSize;
            matrixSize *= (src == dst) ? (srcSize - 1) : insideSupernode[dst].size();
            double weight = (double)tmpEdgeNum / matrixSize;
            tmpDataCost = entropy(weight) * matrixSize;
        }
        else{
            tmpDataCost = sparseEdgeCost * tmpEdgeNum;
        }
        if(src==dst) tmpDataCost/=2;
        return tmpDataCost;
    }

    // Returns minimum total description cost of an edge (src U dst) if they are merged
    protected double dataCostPair(int src, int dst, int n){
        Int2IntOpenHashMap source = superGraph[src];
        Int2IntOpenHashMap dest = superGraph[dst];

        int mergeSize = insideSupernode[src].size() + insideSupernode[dst].size();
        long matrixSize = mergeSize;
        int tmpEdgeNum;

        if(n==src){
            tmpEdgeNum = (source.get(src) & 0x7FFFFFFF) + 2 * (source.get(dst) & 0x7FFFFFFF) + (dest.get(dst) & 0x7FFFFFFF);
            matrixSize *= (matrixSize - 1);
        }
        else{
            tmpEdgeNum = (source.get(n) & 0x7FFFFFFF) + (dest.get(n) & 0x7FFFFFFF);
            matrixSize *= insideSupernode[n].size();
        }
        double weight = (double)tmpEdgeNum / matrixSize;
        if(src==n){
            tmpEdgeNum/=2;
            matrixSize/=2;
        }
        double dense = 2 * log2(numSuperNodes) + log2(maxWeight) + entropy(weight)*matrixSize;
        double sparse = sparseEdgeCost*tmpEdgeNum;

        double tmpMergeCost = (dense<=sparse) ? dense : (-sparse);

        return tmpMergeCost;
    }

    public void removeSupernode(int target){
        numSuperNodes--;
        snList.set(rep[target], snList.getInt(numSuperNodes));
        rep[snList.getInt(rep[target])] = rep[target];
        rep[target] = -1;
        snList.popInt();
    }

    // Merge the best candidate pair
    // Update the dataCost & modelCost of the new merger (src U dst)
    public void merge(int src, int dst){
        int tmpEdgeNum, tmpMax = 0, tmpModelCost, modelTmp;
        int denseNum;
        double tmpMergeCost, tmpDataCost, dataTmp;

        Int2IntOpenHashMap source = superGraph[src];
        Int2IntOpenHashMap dest = superGraph[dst];
        IntOpenHashSet mergeNodeSet = mergeNodeSet(src, dst);

        for(int n: mergeNodeSet){
            denseNum = 0;
            if(n == src){
                tmpEdgeNum = (source.get(src) & 0x7FFFFFFF) + 2 * (source.get(dst) & 0x7FFFFFFF) + (dest.get(dst) & 0x7FFFFFFF);
                if(superGraph[src].get(src) > 0) denseNum++;
                if(superGraph[src].get(dst) > 0) denseNum++;
                if(superGraph[dst].get(dst) > 0) denseNum++;
            }
            else{
                tmpEdgeNum = (source.get(n) & 0x7FFFFFFF) + (dest.get(n) & 0x7FFFFFFF);
            }
            if(tmpMax<tmpEdgeNum) tmpMax = tmpEdgeNum;
            if(n!=src){
                dataTmp = dataCost[n] - (dataCostPair(src,n)+dataCostPair(dst,n));
                modelTmp = modelCost[n];

                if(superGraph[src].get(n) > 0){
                    modelTmp -= 1;
                    denseNum++;
                }

                if(superGraph[dst].get(n) > 0){
                    modelTmp -= 1;
                    denseNum++;
                }
                dataCost[n] = dataTmp;
                modelCost[n] = modelTmp;
            }

            tmpMergeCost = dataCostPair(src, dst, n);

            if(tmpMergeCost>=0){
                numSuperEdges -= (denseNum - 1);
                modelCost[n]++;
            }

            else {
                numSuperEdges -= denseNum;

                tmpEdgeNum |= 0x80000000;
            }
            //Update src tmpEdgeNum
            superGraph[src].put(n, tmpEdgeNum);
            superGraph[n].put(src, tmpEdgeNum);

            //Delete dst from superGraph
            superGraph[n].remove(dst);
        }

        if (tmpMax > maxWeight) maxWeight = tmpMax;

        //Delete dst from superGraph
        superGraph[dst] = null;

        //update inside info & delete supernode dst
        insideSupernode[src].addAll(insideSupernode[dst]);
        insideSupernode[dst] = null;

        //Update tmpDataCost & tmpModelCost
        mergeNodeSet.add(src);
        for(int n: mergeNodeSet){
            if(n == src){
                tmpDataCost = 0;
                tmpModelCost = 0;
                for(int m: superGraph[src].keySet()){
                    tmpDataCost += dataCostPair(src, m);
                    if(superGraph[src].get(m) > 0) tmpModelCost += 1;
                }
                dataCost[src] = tmpDataCost;
                modelCost[src] = tmpModelCost;
            }
            else{
                dataCost[n] += dataCostPair(src,n);
            }
        }
        dataCost[dst] = 0; modelCost[dst] = 0;
        removeSupernode(dst);
    }

    // Caculate MDLsaving of the merger (src & dst)
    public double savingMDL(int src, int dst, double edgeCost){
        double tmpDataCost = dataCost[src] + dataCost[dst];
        double tmpModelCost = edgeCost * (modelCost[src] + modelCost[dst]);
        double denominator, tmpMergeCost = 0, tmpData;
        int dense=0, sparse=0;
        IntOpenHashSet mergeNodeSet = mergeNodeSet(src, dst);

        // Subtract the double counted part
        tmpModelCost -= (superGraph[src].get(dst) > 0) ? edgeCost : 0;
        tmpDataCost -= dataCostPair(src, dst);
        denominator = tmpModelCost + tmpDataCost;

        for(int n: mergeNodeSet){
            tmpData = dataCostPair(src, dst, n);
            if(tmpData<0) tmpMergeCost -= tmpData;
            else tmpMergeCost += tmpData;
        }
        return 1-(tmpMergeCost/denominator);
    }

    public void dropEdges(){
        Int2BooleanOpenHashMap tmp;
        IntOpenHashSet nodeSet = new IntOpenHashSet();
        int check;


        for(int n: snList){
            superGraph[n].int2IntEntrySet().removeIf(e -> (e.getIntValue() < 0));
        }

        try {
            int i = 0;
            while(i < snList.size()){
                int target = snList.getInt(i);
                i++;
                if(superGraph[target].size() == 0) {
                    if (isolatedId == -1) {
                        isolatedId = target;
                    } else {
                        insideSupernode[isolatedId].addAll(insideSupernode[target]);
                        superGraph[target] = null;
                        insideSupernode[target] = null;
                        removeSupernode(target);
                        i--;
                    }
                }
            }
        }catch(Exception e){
            e.printStackTrace();
            System.out.println("Error ");
        }
    }

    public double log2(double x){
        return Math.log(x)/log2;
    }

    // Generate signature matrix
    public int[][] signature(int k){
        final int maxCore = Runtime.getRuntime().availableProcessors();
        ExecutorService executorService = Executors.newFixedThreadPool(maxCore);
        List<Future<int[]>> resultList = new ArrayList<>();
        int[][] ans = new int[numSuperNodes][k+1];


        for (int i=0; i<k; i++)
        {
            MinH calculate = new MinH(numSuperNodes, numNodes, snList, insideSupernode, edges);
            Future<int[]> result = executorService.submit(calculate);
            resultList.add(result);
        }


        executorService.shutdown();
        try {
            executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
        }

        int col = 0;
        for(Future<int[]> future : resultList)
        {
            try
            {
                int j = 0;
                for(int i:future.get()){
                    ans[j++][col] = i;
                }
            }
            catch (InterruptedException | ExecutionException e)
            {
                e.printStackTrace();
            }
            col++;
        }

        for(int i=0;i<numSuperNodes;i++) {
            ans[i][col] = snList.getInt(i);
        }

        int [][] sorted = Arrays.stream(ans).sorted(Arrays::compare).toArray(int[][]::new);

        return sorted;
    }

    public long SSumM(int T, int k) {
        long startTime, endTime;
        double threshold;
        int check=-1;
        SimpleDateFormat format = new SimpleDateFormat( "yyyy-MM-dd HH:mm:ss");
        String formatTime = format.format(System.currentTimeMillis());
        System.out.println("Start Time: " + formatTime);
        startTime = System.currentTimeMillis();
        for(int m = 1; m <= T; m++){
            System.out.println("iter: " + m);
            if(m != T) {
                threshold = 1.0 / (1 + m);
            }else{
                threshold = 0;
            }
            int[][] sorted = signature(k);
            ArrayList<int[]> tmp = divSupernode(sorted, k);

            check = mergeStep(tmp, threshold);
            if(check == 0) break;
        }
        dropEdges();
        if(targetSummarySize >= numNodes * log2(numSuperNodes) + numSuperEdges * (2 * log2(numSuperNodes) + log2(maxWeight))) check=0;
        if(check == 1){
            topNDrop();
        }
        endTime = System.currentTimeMillis();
        return (endTime - startTime);
    }

    // Further sparsification phase
    public void topNDrop(){
        Int2BooleanOpenHashMap tmp;
        IntOpenHashSet nodeSet = new IntOpenHashSet();
        int tmpEdgeNum, src, dst;

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
                    cost[col] = numE * numE / (double) blockSize;
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

    public int mergeStep(ArrayList<int[]> array, double threshold){
        int length;
        int insideCount;
        double maxSaving, saving, summarySize;
        double superedgeCost = 2 * log2(numSuperNodes) + log2(maxWeight);
        int sampleSize;
        for (int[] tmp : array) {
            length = tmp.length;

            if (length < 2) continue;
            insideCount = 0;
            int l = length - 1, until = 0;

            sampleSize = length;
            while(l>0){l >>= 1; until++;}

            int m = until - 1;
            int n =  sampleSize - 1;
            until = 0;
            sampleSize = 0;
            while(m>0){m >>= 1; until++;}
            while(n>0){n >>= 1; sampleSize++;}

            while (insideCount < until) {
                cntFlag += 1;
                long uniqueValue = cntFlag;
                if (length < 2) break;
                int bestSrc = -1, bestDst = -1;
                maxSaving = -1;
                for(int i = 0; i < length; i++) {
                    int rn = rnd.nextInt(length * (length - 1));
                    int src = rn / (length - 1);
                    int dst = rn % (length - 1);
                    if(src <= dst) dst++;
                    if (dup[src][dst] == uniqueValue) {
                        continue;
                    }
                    dup[src][dst] = uniqueValue;

                    saving = savingMDL(tmp[src], tmp[dst], superedgeCost);
                    if (maxSaving < saving) {
                        maxSaving = saving;
                        bestSrc = src;
                        bestDst = dst;
                    }
                }
                if (maxSaving < threshold) {
                    insideCount++;
                    continue;
                }
                insideCount = 0;
                merge(tmp[bestSrc], tmp[bestDst]);
                length--;
                tmp[bestDst] = tmp[length];
                superedgeCost = 2 * log2(numSuperNodes) + log2(maxWeight);
                summarySize = numNodes * log2(numSuperNodes) + numSuperEdges * (2 * log2(numSuperNodes) + log2(maxWeight));
                if(summarySize <= targetSummarySize) return 0;
            }
        }
        return 1;
    }

    // Divide signature matrix into supernode groups (spec. 500)
    public ArrayList<int[]> divSupernode(int[][] array, int k){
        ArrayList<int[]> tmp;


        int[] supernode = new int[numSuperNodes];


        for(int i=0;i<numSuperNodes;i++){
            supernode[i] = array[i][0];
        }

        tmp = listN(supernode, 0);

        return recursiveDiv(array, tmp, 0, k);
    }

    public ArrayList<int[]> recursiveDiv(int[][] array, ArrayList<int[]> tmp, int col, int k){
        ArrayList<int[]> ans = new ArrayList<>();

        for(int[] x:tmp){
            int first = x[0];
            int last = x[1];
            int length = last - first + 1;
            if((length>divThreshold)&&(col<k)){
                col++;
                int [] tempp = column(array, col, first, last);

                ArrayList<int[]> tempo = listN(tempp, first);

                ans.addAll(recursiveDiv(array, tempo ,col, k));
                col--;
            }

            // Groups are not smaller than divThreshold
            if((length>divThreshold)&&(col==k)){
                // Divide each group by using first & Last
                int quot = length / divThreshold;
                int tmpLast = first + divThreshold -1;
                quot--;
                while(tmpLast<last){
                    ans.add(column(array, k, first, tmpLast));
                    first += divThreshold;
                    if(quot!=0){
                        tmpLast += divThreshold;
                        quot--;
                        continue;
                    }
                    tmpLast = last;
                }
                // Left over
                ans.add(column(array, k, first, tmpLast));
            }
            if(length<=divThreshold) {
                ans.add(column(array, k, first, last));
            }
        }
        return ans;
    }

    public int[] column(int[][] array, int col, int first, int last){
        int[] ans = new int[last-first+1];
        int index = 0;
        for(int i = first; i<=last; i++){
            ans[index++] = array[i][col];
        }
        return ans;
    }

    public ArrayList<int[]> listN(int[] array, int first){
        ArrayList<int[]> ans = new ArrayList<>();
        int tmp = array[0];
        int last = first-1;

        for(int i=0;i<array.length;i++){
            if(tmp==array[i]) last++;
            else{
                int[] model = new int[2];
                model[0] = first;
                model[1] = last;
                ans.add(model);
                tmp=array[i];
                first = last++ +1;
            }
        }
        int[] model = new int[2];
        model[0] = first;
        model[1] = last;
        ans.add(model);
        return ans;
    }

    public double entropy(double p) {
        double entropy;
        if ((p == 0.0) || (p == 1.0)) {
            entropy = 0;
        } else {
            entropy = -(p * log2(p) + (1 - p) * log2(1 - p));
        }
        return entropy;
    }

    public IntOpenHashSet mergeNodeSet(int src, int dst){
        IntOpenHashSet returnSet = new IntOpenHashSet();
        returnSet.addAll(superGraph[src].keySet());
        returnSet.addAll(superGraph[dst].keySet());
        if(returnSet.contains(dst)){
            returnSet.add(src);
            returnSet.remove(dst);
        }
        return returnSet;
    }

    public static void main(String[] args) {
        String filename = args[0];
        double fracK = Double.parseDouble(args[1]);
        int re = Integer.parseInt(args[2]);
        final SSumM graph;
        fracK = Math.ceil(fracK*100)/100.0;

        String dataPath = filename;
        System.out.println("---------------------------------------------------");
        if(re == 1) {
            graph = new ReOne(dataPath);
        }else{
            graph = new SSumM(dataPath);
        }
        graph.inputGraph();
        graph.targetSummarySize = fracK * 2 * graph.numEdges * graph.log2(graph.numNodes);
        double elapsTime = graph.SSumM(20, 10);
        graph.summarySave(filename, fracK);
        double error = graph.norm();
        double originalSize = 2 * graph.numEdges * graph.log2(graph.numNodes);
        double summarySize = graph.numNodes * graph.log2(graph.numSuperNodes) + graph.numSuperEdges * (2 * graph.log2(graph.numSuperNodes) + graph.log2(graph.maxWeight));
        double frac = summarySize/originalSize;
        System.out.println("---------------------------------------------------");
        System.out.println("Elapsed Time: " + elapsTime+" ms");
        System.out.println("Original size: " + String.format("%.2f", originalSize) + " bits");
        System.out.println("Summary size: " + String.format("%.2f", summarySize) + " bits (" + String.format("%.6f", 100*summarySize/originalSize) + "%)");
        System.out.println("L"+re+"error: " +  String.format("%.2e", error));
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
            int edgeCnt = superGraph[u].getOrDefault(v, 0);
            long sz = insideSupernode[u].size();
            sz *= ((u == v) ? (insideSupernode[u].size() - 1) : insideSupernode[v].size());
            double w = edgeCnt / (double)sz;
            err += ((1-w) * (1-w) - w * w);
        }

        err *= 2;
        for(int v: snList){
            for(Int2IntMap.Entry u: superGraph[v].int2IntEntrySet()){
                long sz = insideSupernode[u.getIntKey()].size();
                sz *= ((u.getIntKey() == v) ? (insideSupernode[u.getIntKey()].size() - 1) : insideSupernode[v].size());
                int edgeCnt = u.getIntValue();
                double w = edgeCnt / (double)sz;
                err += (w * w * sz);
            }
        }
        err = Math.sqrt(err);
        err /= numNodes;
        err /= (numNodes - 1);
        return err;
    }

    public void summarySave(String filename, double fracK) {
        Int2IntOpenHashMap idx2Node = new Int2IntOpenHashMap();
        for (Int2IntMap.Entry nti: node2Idx.int2IntEntrySet()){
            idx2Node.put(nti.getIntValue(), nti.getIntKey());
        }
        String[] fileNames = filename.split("/");
        String file = fileNames[fileNames.length -1].split("\\.")[0];
        String algName = "";

        filename = "summary_" + file + ".txt";
        File f = new File("output/" + filename);
        try {
            FileWriter fw = new FileWriter(f);
            fw.write("<Subnode of each supernode>");
            fw.write(System.getProperty( "line.separator" ));
            for (int sup_v : snList) {
                fw.write(String.format("%d", sup_v));
                for (int sub_v : insideSupernode[sup_v]) {
                    fw.write("\t"+String.format("%d", idx2Node.get(sub_v)));
                }
                fw.write(System.getProperty( "line.separator" ));
            }
            fw.write("<Superedge info>");
            fw.write(System.getProperty( "line.separator" ));
            for(int sup_v: snList){
                for(int neighbor_v : superGraph[sup_v].keySet()){
                    if(sup_v>=neighbor_v){
                        if(sup_v == neighbor_v) fw.write(String.format("%d", sup_v) +"\t"+ String.format("%d", neighbor_v) +"\t"+ String.format("%d", superGraph[sup_v].get(neighbor_v)/2));
                        else{
                            fw.write(String.format("%d", sup_v) +"\t"+ String.format("%d", neighbor_v) +"\t"+ String.format("%d", superGraph[sup_v].get(neighbor_v)));
                        }
                        fw.write(System.getProperty( "line.separator" ));
                    }
                }
            }
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
