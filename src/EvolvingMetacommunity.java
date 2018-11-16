

import java.util.Random;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;


enum SelAction 
{
	REPR, MORT
}

public class EvolvingMetacommunity {

	public static void main(String[] args) throws IOException {
		
		Comm comm = new Comm();
		Evol evol = new Evol();
		Run run = new Run();
		if (args.length > 0)
			Reader.readInput(args[0], comm, evol, run);
		comm.init();
		evol.init();

		int processors = Runtime.getRuntime().availableProcessors();
		System.out.println(Integer.toString(processors) + " processor"
				+ (processors != 1 ? "s are " : " is ")
				+ "available");
		Run.processors = processors; 

		try (PrintWriter streamOut = new PrintWriter(new FileWriter(run.fileName))) {

			long startTime = System.currentTimeMillis();

			streamOut.println("gridsize;nbr_patches;p_ext;m;ro;nbr_envs;sigma_env;microsites;K;d;b;loci;sigma_ftp;mu;w_e;sex;"
					+ "run;time;patch;e_p;" 
					+ "species;N;genotype;genotype_var;fenotype;fenotype_var;fitness;fitness_var");

			for (int r = 0; r < run.runs; r++) {

				System.out.format("run = %d%n", (r+1));

				Init init = new Init(comm, evol);
				Sites sites = new Sites(comm, evol, init);

				System.out.format("  time = %d; metacommunity N = %d%n", 0, sites.metacommunitySize());

				for (int t = 0; t < run.timeSteps; t++) {

// debug					
//			        System.out.println("  " + java.util.Arrays.toString(sites.alive));
//
			        
					sites.patchExtinction();

					sites.mortality();

// debug					
//			    	System.out.println("  mortality");
//			        System.out.println("  ID : " + java.util.Arrays.toString(sites.ID));
//			        System.out.println("  alive : " + java.util.Arrays.toString(sites.alive));
//			        System.out.println("  posAdults : " + java.util.Arrays.toString(sites.posAdults));
//			        System.out.println("  nbrAdults : " + java.util.Arrays.toString(sites.nbrAdults));
//			        System.out.println("  endPosAdults : " + sites.endPosAdults);
//			        System.out.println("  cumsumAdults : " + java.util.Arrays.toString(sites.cumsumAdults));
//			        System.out.println("  posEmpty : " + java.util.Arrays.toString(sites.posEmpty));
//			        System.out.println("  nbrEmpty : " + java.util.Arrays.toString(sites.nbrEmpty));
//			        System.out.println("  cumsumEmpty : " + java.util.Arrays.toString(sites.cumsumEmpty));
//

			        sites.reproduction();

// debug					
//			    	System.out.println("  reproduction");
//			        System.out.println("  " + java.util.Arrays.toString(sites.alive));
//			    	System.out.println("");
//

					if (t == 0 || ((t+1) % run.printSteps) == 0) 
					{
						System.out.format("  time = %d; metacommunity N = %d%n", (t+1), sites.metacommunitySize());
					}
					if (t == 0 || ((t+1) % run.saveSteps) == 0) 
					{
						for (int p = 0; p < comm.nbrPatches; p++) 
							for (int s = 0; s < comm.nbrSpec; s++ )
							{
								streamOut.format("%d;%d;%f;%f;%f;%d;%f;%d;%d;%f;%f;%d;%f;%f;%f;%s",
										comm.gridSize,comm.nbrPatches,comm.pExt,comm.dispRate,comm.ro,comm.nbrEnv,comm.envStd,comm.microsites,comm.K,comm.d,comm.b,evol.loci,evol.ftpStd,evol.mutationRate,evol.wE,comm.sex);
								streamOut.format(";%d;%d;%d;%f",
										r+1,t+1,p+1,sites.patchEnv[p]);
								streamOut.format(";%d;%d;%f;%f;%f;%f;%f;%f%n",
										s+1,sites.popSize(p, s),sites.specGenotype(p,s),sites.specGenotypeVar(p,s),sites.specFenotype(p,s),sites.specFenotypeVar(p,s),sites.specFitness(p,s),sites.specFitnessVar(p,s));
							}
					}
				}
			}

			long endTime = System.currentTimeMillis();
			System.out.println("EvolMetac took " + (endTime - startTime) +
					" milliseconds.");

			streamOut.close();
		}
	}
}


class Sites 
{
	static final ForkJoinPool pool = new ForkJoinPool();
    Comm comm;
    Evol evol;
    int totSites;
    
    int[] patch;
    double[] resource;
    boolean[] alive;
    int[] ID;
    double[] fenotp;
    double[] fitness;
    int[][] genotp;

    int[] posAdults;
    int[][] posAdultsPatch;
    int[] nbrAdults;
    int[] cumsumAdults;
    int endPosAdults;
    int[] endPosAdultsPatch;
    int[] posEmpty;
    int[] nbrEmpty;
    int[] cumsumEmpty;
    double[] adultsContr;

    double[] patchEnv;
    
    public Sites (Comm cmm, Evol evl, Init init) 
    {
        comm = cmm;
        evol = evl;
        totSites = comm.nbrPatches * comm.microsites;
        
        comm.calcDispNeighbours();
        
        patch = new int [totSites];
        resource = new double [totSites];
        alive = new boolean [totSites];
        ID = new int [totSites];
        fenotp = new double [totSites];
        fitness = new double [totSites];
        genotp = new int [totSites][2*evol.loci];

        posAdults = new int[totSites];
        posAdultsPatch = new int[comm.nbrPatches][comm.microsites];
        nbrAdults = new int[comm.nbrPatches];
        cumsumAdults = new int[comm.nbrPatches];
        endPosAdultsPatch = new int[comm.nbrPatches];
        posEmpty = new int[totSites];
        nbrEmpty = new int[comm.nbrPatches];
        cumsumEmpty = new int[comm.nbrPatches];
        adultsContr = new double[totSites];
        
        patchEnv = new double[comm.nbrPatches];
        
        double indGtp = 0.5;
		int maxOnesMother = 0;
		int maxOnesFather = 0;

		java.util.Arrays.fill(alive,false);
 
        for (int p = 0; p < comm.nbrPatches; p++) {
        	patchEnv[p] = init.resource[p];
            for (int m = 0; m < comm.microsites; m++) {
                patch[(p * comm.microsites) + m] = p;
                resource[(p * comm.microsites) + m] = patchEnv[p] + (Auxils.random.nextGaussian() * comm.envStd);
            }
            int[] posInds = Auxils.arraySample(init.N[p], Auxils.enumArray(p * comm.microsites, ((p + 1) * comm.microsites) - 1));

// debug					
//	        System.out.println("patch  " + p);
//	        System.out.println("  init resource  " + init.resource[p]);
//	        System.out.println("  patchEnv  " + patchEnv[p]);
//	        System.out.println("  init genotp  " + init.genotp[p]);
//	    	System.out.println("");
//

	    	for (int m : posInds) {
            	alive[m] = true;
            	ID[m] = init.ID[p];
            	fitness[m] = 1;

            	indGtp = init.genotp[p];
            	
            	maxOnesMother = (int) (indGtp*evol.loci);
            	maxOnesMother = maxOnesMother > evol.loci ? evol.loci : maxOnesMother; 
            	maxOnesMother = maxOnesMother < 0 ? 0 : maxOnesMother; 
            	maxOnesFather = (int) (indGtp*evol.loci + 0.5);
            	maxOnesFather = maxOnesFather > evol.loci ? evol.loci : maxOnesFather; 
            	maxOnesFather = maxOnesFather < 0 ? 0 : maxOnesFather; 
            	for (int l = 0; l < maxOnesMother; l++) {
            		genotp[m][evol.allMother[l]] = 1;
            	}
            	for (int l = 0; l < maxOnesFather; l++) {
            		genotp[m][evol.allFather[l]] = 1;
            	}
            	
// debug					
//		    	System.out.println("    ind  " + m);
//		        System.out.println("      " + java.util.Arrays.toString(genotp[m]));
//		    	System.out.println("");
//

            	fenotp[m] = Auxils.arrayMean(Auxils.arrayElements(genotp[m],evol.allTrait)) + (Auxils.random.nextGaussian() * evol.ftpStd);
            	fitness[m] = Math.exp(-(Math.pow(fenotp[m]-resource[m],2))/evol.divK);
            }
        }
    }
    
    void patchExtinction()
    {
        for (int p = 0; p < comm.nbrPatches; p++)
            if (Auxils.random.nextDouble() <= comm.pExt)
                for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
                    alive[i] = false;
    }
            
    void mortality()
    {
        int countDead = 0;
        endPosAdults = 0;
        java.util.Arrays.fill(nbrAdults, 0);
        java.util.Arrays.fill(nbrEmpty, 0);
        for (int i = 0; i < alive.length; i++) {
        	if (alive[i]) {
	        	switch (evol.selAction) {
	        		case REPR:
	        			alive[i] = Auxils.random.nextDouble() > comm.d;
	        			break;
	        		case MORT: 
	        			alive[i] = Auxils.random.nextDouble() < (1-comm.d)*fitness[i];
	        			break;
	        	}
        	}
            if (alive[i]) {
                posAdults[endPosAdults++] = i;
                nbrAdults[patch[i]]++;
            }
            else {
                posEmpty[countDead++] = i;
                nbrEmpty[patch[i]]++;
            }
        }
        System.arraycopy(nbrEmpty, 0, cumsumEmpty, 0, nbrEmpty.length);
        Auxils.arrayCumSum(cumsumEmpty);
        System.arraycopy(nbrAdults, 0, cumsumAdults, 0, nbrAdults.length);
        Auxils.arrayCumSum(cumsumAdults);
    }

    void reproduction()
    {   
    	contributionAdults();
        ForkReproduction.setThreshold(20);
        ForkReproduction fr = new ForkReproduction(this, 0, comm.nbrPatches);
        pool.invoke(fr);
    }    
    
    void contributionAdults()
    {
    	ForkContribution.setThreshold(10000);
    	ForkContribution ff = new ForkContribution(this, 0, endPosAdults);
        pool.invoke(ff);
    }
    
    int metacommunitySize()
    {
    	int sum = Auxils.arraySum(alive);
        return sum;
    }

    int communitySize(int p)
    {
        int sum = Auxils.arraySum(java.util.Arrays.copyOfRange(alive, p*comm.microsites, p*comm.microsites + comm.microsites));
        return sum;
    }

    int popSize(int p, int s)
    {
        int sum = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i] && ID[i] == s)
        		sum++;
        return sum;
    }

    double meanEnv(int p)
    {
        double mean = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
    		mean += resource[i];
            mean /= comm.microsites;
            return mean;
    }

    double specGenotype(int p, int s)
    {
        double mean = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i] && ID[i] == s)
        		mean += Auxils.arrayMean(Auxils.arrayElements(genotp[i],evol.allTrait));
        mean /= popSize(p,s);
        return mean;
    }

    double specGenotypeVar(int p, int s)
    {
        double mean = specGenotype(p, s);
        double var = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i] && ID[i] == s)
        		var += Math.pow(mean - Auxils.arrayMean(Auxils.arrayElements(genotp[i],evol.allTrait)), 2);
        var /= popSize(p,s);
        return var;
    }

    double specFenotype(int p, int s)
    {
        double mean = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i] && ID[i] == s)
        		mean += fenotp[i];
        mean /= popSize(p,s);
        return mean;
    }

    double specFenotypeVar(int p, int s)
    {
        double mean = specFenotype(p, s);
        double var = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i] && ID[i] == s)
        		var += Math.pow(mean - fenotp[i], 2);
        var /= popSize(p,s);
        return var;
    }

    double specFitness(int p, int s)
    {
        double mean = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i] && ID[i] == s)
        		mean += fitness[i];
        mean /= popSize(p,s);
        return mean;
    }

    double specFitnessVar(int p, int s)
    {
        double mean = specFitness(p, s);
        double var = 0;
        for (int i = p*comm.microsites; i < (p+1)*comm.microsites; i++)
        	if (alive[i] && ID[i] == s)
        		var += Math.pow(mean - fitness[i], 2);
        var /= popSize(p,s);
        return var;
    }
}


class ForkContribution extends RecursiveAction 
{   
	private static final long serialVersionUID = 1L;
	protected static int sThreshLim = 100;
    protected static int sThreshold = sThreshLim;
	Sites sites;
    int iStart, iNbr;
    
    public ForkContribution(Sites s, int is, int in)
    {
    	sites = s;
    	iStart = is;
    	iNbr = in;
    }

    public static void setThreshold(int tr)
    {
    	if (tr > sThreshLim)
    		sThreshold = tr;
    }

    protected void computeDirectly() 
    {
    	for (int i = iStart; i < iStart + iNbr; i++) 
    	{
			if (sites.nbrAdults[sites.patch[sites.posAdults[i]]] < sites.comm.K)
				sites.adultsContr[i] = (1 - (double)sites.nbrAdults[sites.patch[sites.posAdults[i]]]/sites.comm.K);
			else
				sites.adultsContr[i] = 0;
			if (sites.evol.selAction == SelAction.REPR)
				sites.adultsContr[i] *= sites.fitness[sites.posAdults[i]];
    	}

// debug
//    	System.out.println("  contribution");
//    	System.out.println("  K : " + sites.comm.K);
//        System.out.println("  nbrAdults : " + java.util.Arrays.toString(sites.nbrAdults));
//        System.out.println("  contribution : " + java.util.Arrays.toString(sites.adultsContr));
//
        
    }
    
    protected void compute() 
    {
    	if (iNbr < sThreshold) {
    		computeDirectly();
    		return;
    	}

    	int split = iNbr / 2;

    	invokeAll(new ForkContribution(sites, iStart, split),
    			new ForkContribution(sites, iStart + split, iNbr - split));
    }
}


class ForkReproduction extends RecursiveAction 
{   
	private static final long serialVersionUID = 1L;
	protected static int sThreshLim = 2;
    protected static int sThreshold = sThreshLim;
	Sites sites;
    int pStart, pNbr;

    int nbrSettled, startPosSample, rangeSample, posFather, nbrFathers, patchMother;
    double nbrOffspring;
    int[] posOffspring, posMothers, sampleM;
    double[] adultsProb;
    
    public ForkReproduction(Sites s, int ps, int pn)
    {
    	sites = s;
    	pStart = ps;
    	pNbr = pn;
        adultsProb = new double[sites.totSites];
    }
   
    public static void setThreshold(int tr)
    {
    	if (tr > sThreshLim)
    		sThreshold = tr;
    }

    protected void computeDirectly() 
    {

// debug
//		System.out.println("  reproduction");
//		
		
    	for (int p = pStart; p < pStart + pNbr; p++) {
    		probParents(p);
    		nbrOffspring = Auxils.arraySum(adultsProb, sites.endPosAdults)*sites.comm.b;
	    	if (nbrOffspring < 1)
    			nbrOffspring = (Auxils.random.nextDouble() <= nbrOffspring) ? 1 : 0;

// debug
//	    	System.out.println("  patch : " + p);
//	    	System.out.println("  probParents : " + java.util.Arrays.toString(adultsProb));
//	    	System.out.println("  nbrOffspring : " + nbrOffspring);
//	    	    		
	    		    	
    		nbrSettled = Math.min((int) (nbrOffspring + 0.5), sites.cumsumEmpty[p] - ((p == 0) ? 0 : sites.cumsumEmpty[p-1]));
    		if (nbrSettled > 0) {
    			posOffspring = Auxils.arraySample(nbrSettled, java.util.Arrays.copyOfRange(sites.posEmpty, (p == 0) ? 0 : sites.cumsumEmpty[p-1], sites.cumsumEmpty[p]));
    			sampleM = Auxils.arraySampleProb(nbrSettled, Auxils.enumArray(0, sites.endPosAdults-1), java.util.Arrays.copyOf(adultsProb, sites.endPosAdults));
    			posMothers = Auxils.arrayElements(sites.posAdults, sampleM);

// debug
//    	    	System.out.println("  nbrSettled : " + nbrSettled);
//    	        System.out.println("  posOffspring : " + java.util.Arrays.toString(posOffspring));
//    	        System.out.println("  sampleM : " + java.util.Arrays.toString(sampleM));
//    	        System.out.println("  posMothers : " + java.util.Arrays.toString(posMothers));
//

    			settle(p);
    		}
    	}
    }
    
    protected void compute() 
    {
        if (pNbr < sThreshold) {
            computeDirectly();
            return;
        }
 
        int split = pNbr / 2;
 
        invokeAll(new ForkReproduction(sites, pStart, split),
                new ForkReproduction(sites, pStart + split, pNbr - split));
    }
    
    void probParents(int p)
    {
    	for (int i = 0; i < sites.endPosAdults; i++)
    		adultsProb[i] = sites.adultsContr[i]*sites.comm.dispNeighbours[p][sites.patch[sites.posAdults[i]]];
    }
    
    void settle(int p)
    {
        int posFather, nbrFathers, patchMother;
        int[] posFathers;
        double[] FathersProb;

// debug
//    	System.out.println("  settle");
//    	System.out.println("  comm.sex : " + sites.comm.sex);
//
        
    	for (int i = 0; i < nbrSettled; i++) {
    		sites.alive[posOffspring[i]] = true;
    		sites.ID[posOffspring[i]] = sites.ID[posMothers[i]];
    		sites.fitness[posOffspring[i]] = 1;

// debug
//    		System.out.println("  offspr : " + i);
//    		System.out.println("  posOffspr : " + posOffspring[i]);
//    		System.out.println("  posMother : " + posMothers[i]);
//    		System.out.println("  patch mother : " + sites.patch[posMothers[i]]);
//    		System.out.println("  ID mother : " + sites.ID[posMothers[i]]);
//    		System.out.println("  ID offspring : " + sites.ID[posOffspring[i]]);
//
    		    		        
    		if (sites.comm.sex) {

// debug
//    			System.out.println("  in sex");
//
    			    		        
    			patchMother = sites.patch[posMothers[i]];
    			startPosSample = (patchMother == 0) ? 0 : sites.cumsumAdults[patchMother-1];
    			rangeSample = sites.cumsumAdults[patchMother] - startPosSample; 
    			posFathers = new int[rangeSample];
    			FathersProb = new double[rangeSample];
    			nbrFathers = 0;
    			for (int f = startPosSample; f < (startPosSample + rangeSample); f++)
    				if (sites.ID[sites.posAdults[f]] == sites.ID[posMothers[i]]) {
    					posFathers[nbrFathers] = sites.posAdults[f];
    					FathersProb[nbrFathers++] = adultsProb[f];
    				}
    			posFather = posFathers[Auxils.randIntProb(nbrFathers, FathersProb)];

// debug
//    			System.out.println("  patch mother : " + patchMother);
//    			System.out.println("  patch father : " + sites.patch[posFather]);
//    			System.out.println("  posMother : " + posMothers[i]);
//    			System.out.println("  posFather : " + posFather);
//    			System.out.println("  ID mother : " + sites.ID[posMothers[i]]);
//    			System.out.println("  ID father : " + sites.ID[posFather]);
//
    			    			    		        
    			inherit(posOffspring[i], posMothers[i], posFather);
    		}
    		else
    			inherit(posOffspring[i], posMothers[i]);

    		sites.fenotp[posOffspring[i]] = Auxils.arrayMean(Auxils.arrayElements(sites.genotp[posOffspring[i]],sites.evol.allTrait)) + (Auxils.random.nextGaussian() * sites.evol.ftpStd);
    		sites.fitness[posOffspring[i]] = Math.exp(-(Math.pow(sites.fenotp[posOffspring[i]]-sites.resource[posOffspring[i]],2))/sites.evol.divK);
    	}
    }
    
    void inherit(int posOffspring, int posParent)
    {
    	for (int l : sites.evol.allTrait) {
    		if (Auxils.random.nextDouble() <= sites.evol.mutationRate)
    			sites.genotp[posOffspring][l] = (sites.genotp[posParent][l] == 0) ? 1 : 0;
    		else
    			sites.genotp[posOffspring][l] = sites.genotp[posParent][l];
    	}
    }

    void inherit(int posOffspring, int posMother, int posFather)
    {
    	for (int l = 0; l < sites.evol.loci; l++) {
    		if (Auxils.random.nextBoolean())
    			sites.genotp[posOffspring][sites.evol.allMother[l]] = sites.genotp[posMother][sites.evol.allMother[l]];
    		else
    			sites.genotp[posOffspring][sites.evol.allMother[l]] = sites.genotp[posMother][sites.evol.allFather[l]];
    		if (Auxils.random.nextDouble() <= sites.evol.mutationRate)
    			sites.genotp[posOffspring][sites.evol.allMother[l]] = (sites.genotp[posOffspring][sites.evol.allMother[l]] == 0) ? 1 : 0;
    		if (Auxils.random.nextBoolean())
    			sites.genotp[posOffspring][sites.evol.allFather[l]] = sites.genotp[posFather][sites.evol.allMother[l]];
    		else
    			sites.genotp[posOffspring][sites.evol.allFather[l]] = sites.genotp[posFather][sites.evol.allFather[l]];
    		if (Auxils.random.nextDouble() <= sites.evol.mutationRate)
    			sites.genotp[posOffspring][sites.evol.allFather[l]] = (sites.genotp[posOffspring][sites.evol.allFather[l]] == 0) ? 1 : 0;
    	}
    }
}


class Comm 
{
    int nbrEnv = 2;
    int nbrSpec = nbrEnv;
    boolean propExpl = false;
    double propFirst = 1.0/nbrEnv;
    double minResource = 0.4;
    double stepResource = 0.2;
    double envStd = 0.02;
    int microsites = 600;
    int K = 400;
    double d = 0.1;
    double b = 2;
    int gridSize = 2;
    int nbrPatches = gridSize*gridSize; 
    double pExt = 1e-3;
    double dispRate = 0.01;
    double ro = 2;
    boolean sex = true;
    
    double[][] neighbours = new double[nbrPatches][nbrPatches];
    double[][] dispNeighbours = new double[nbrPatches][nbrPatches];

    void init()
    {
    	if (!propExpl)
    		propFirst = 1.0/nbrEnv;
        nbrPatches = gridSize*gridSize;
        neighbours = new double[nbrPatches][nbrPatches];
        dispNeighbours = new double[nbrPatches][nbrPatches];
        calcDistNeighbours();
    }
    
    void calcDistNeighbours() 
    {
        for (int i = 0; i < gridSize; i++)
            for (int j = 0; j < gridSize; j++)
                for (int i2 = 0; i2 < gridSize; i2++)
                    for (int j2 = 0; j2 < gridSize; j2++) {
                       double dist = Math.sqrt(Math.pow(Math.min(Math.abs(i - i2),gridSize-Math.abs(i - i2)),2) + Math.pow(Math.min(Math.abs(j - j2),gridSize-Math.abs(j - j2)),2));
                       neighbours[j*gridSize+i][j2*gridSize+i2] = dist;
                    }
    }
    
    void calcDispNeighbours()
    {
        for (int i = 0; i < nbrPatches; i++) {
            for (int j = 0; j < nbrPatches; j++)
                dispNeighbours[i][j] = (i == j) ? 0 : (ro*Math.exp(-ro*neighbours[i][j]));
            double iSum = Auxils.arraySum(dispNeighbours[i]);
            for (int j = 0; j < nbrPatches; j++)
                dispNeighbours[i][j] = (i == j) ? (1 - dispRate) : (dispRate*dispNeighbours[i][j]/iSum);
        }
    }
}


class Evol 
{
    double wE = 0.2; 
    double divK = 1; 
    int loci = 10;
    double mutationRate = 1e-4;
    double ftpStd = 0.02;
    
    SelAction selAction = SelAction.REPR;

    int [] allMother;
    int [] allFather;
    int [] allTrait;
            
    void init() 
    {
        allMother  = new int [loci];
        allFather  = new int [loci];
        allTrait = new int[2*loci];

        divK = 2*Math.pow(wE, 2);
        
        for (int l = 0; l < loci; l++) {
            allMother[l] = l;
            allFather[l] = l + loci;
        }
        allTrait = Auxils.arrayConcat(allMother, allFather);
    }
}


class Run 
{
	static int processors = 1;
    int runs = 1;
    int timeSteps = 10000;
    int printSteps = 100;
    int saveSteps = 1000;
    String fileName = "output_evolvingMetacommunity.csv";
}


class Init 
{
    int[] theShuffle;
	double[] resource;
    int[] ID;
    int[] N;
    int[] specPatch;
    double[] genotp;
    
    public Init(Comm comm, Evol evol) 
    {
		int sp, env;
    	theShuffle = Auxils.enumArray(0,comm.nbrPatches-1);
        Auxils.arrayShuffle(theShuffle);
    	resource = new double [comm.nbrPatches];
        ID = new int [comm.nbrPatches];
        N = new int [comm.nbrPatches];
		genotp = new double [comm.nbrPatches];

		java.util.Arrays.fill(N, comm.K);

		//		
		sp = 0;
        env = 0;
   		for (int p = 0; p < (int)(comm.propFirst*comm.nbrPatches); p++) 
   		{
   			ID[p] = sp;
        	resource[p] = comm.minResource;
   		}
		sp = 1;
        env = 1;
   		for (int p = (int)(comm.propFirst*comm.nbrPatches); p < comm.nbrPatches; p++) 
   		{
   			if (sp == comm.nbrSpec)
   				sp = 1;
        	if (env == comm.nbrEnv) 
        		env = 1;
   			ID[p] = sp++;
        	resource[p] = (comm.minResource + comm.stepResource*env++);
   		}
		ID = Auxils.arrayElements(ID, theShuffle);
        resource = Auxils.arrayElements(resource, theShuffle);
        genotp = java.util.Arrays.copyOf(resource, resource.length);
    }
}


class Reader
{
	static void readInput(String fileName, Comm comm, Evol evol, Run run) throws IOException
	{
		try (BufferedReader input = new BufferedReader(new FileReader(fileName))) {
			String line;
			String[] words;
			while ((line = input.readLine()) != null) {
				words = line.trim().split("\\s+");
				switch (words[0]) {
					case "COMM": inputComm(words, comm);
					break;
					case "EVOL": inputEvol(words, evol);
					break;
					case "RUN": inputRun(words, run);
					break;
				}
			}

			input.close();
		}
	}

	static void inputComm(String[] in, Comm comm)
	{
		switch (in[1]) {
			case "NBRENV": 
				comm.nbrEnv = Integer.parseInt(in[2]);
				comm.nbrSpec = comm.nbrEnv;
			break;
			case "MINENV": comm.minResource = Double.parseDouble(in[2]);
			break;
			case "STEPENV": comm.stepResource = Double.parseDouble(in[2]);
			break;
			case "PROPFIRST": 
				comm.propFirst = Double.parseDouble(in[2]);
				comm.propExpl = true;
			break;
			case "MICROSITES": comm.microsites = Integer.parseInt(in[2]);
			break;
			case "K": comm.K = Integer.parseInt(in[2]);
			break;
			case "D": comm.d = Double.parseDouble(in[2]);
			break;
			case "B": comm.b = Double.parseDouble(in[2]);
			break;
			case "SEX": 
				switch (in[2]) {
					case "YES": comm.sex = true;
					break;
					case "NO": comm.sex = false;
					break;
				}
			break;
			case "ENVSTD": comm.envStd = Double.parseDouble(in[2]);
			break;
			case "GRIDSIZE": comm.gridSize = Integer.parseInt(in[2]);
			break;
			case "PEXT": comm.pExt = Double.parseDouble(in[2]);
			break;
			case "M": comm.dispRate = Double.parseDouble(in[2]);
			break;
			case "RO": comm.ro = Double.parseDouble(in[2]);
			break;
		}
	}

	static void inputEvol(String[] in, Evol evol)
	{
		switch (in[1]) {
			case "WE": evol.wE = Double.parseDouble(in[2]);
			break;
			case "LOCI": evol.loci = Integer.parseInt(in[2]);
			break;
			case "MU": evol.mutationRate = Double.parseDouble(in[2]);
			break;
			case "FTPSTD": evol.ftpStd = Double.parseDouble(in[2]);
			break;
			case "SELACTION": 
				switch (in[2]) {
					case "REPR": evol.selAction = SelAction.REPR;
					break;
					case "MORT": evol.selAction = SelAction.MORT;
					break;
				}
			break;
		}        
	}

	static void inputRun(String[] in, Run run)
	{
		switch (in[1]) {
			case "RUNS": run.runs = Integer.parseInt(in[2]);
			break;
			case "TIMESTEPS": run.timeSteps = Integer.parseInt(in[2]);
			break;
			case "PRINTSTEPS": run.printSteps = Integer.parseInt(in[2]);
			break;
			case "SAVESTEPS": run.saveSteps = Integer.parseInt(in[2]);
			break;
			case "OUTPUT": run.fileName = in[2];
			break;
		}        
	}
}


class Auxils 
{
    static Random random = new Random();

    static void arrayShuffle(int[] array)
    {
        int index, temp;
        for (int i = array.length - 1; i > 0; i--)
        {
            index = random.nextInt(i + 1);
            temp = array[index];
            array[index] = array[i];
            array[i] = temp;
        }
    }

    static void arrayShuffle(double[] array)
    {
        int index;
        double temp;
        for (int i = array.length - 1; i > 0; i--)
        {
            index = random.nextInt(i + 1);
            temp = array[index];
            array[index] = array[i];
            array[i] = temp;
        }
    }
    
    static int[] arraySample(int n, int[] array)
    {
        int[] tempArr = array.clone();
        arrayShuffle(tempArr);
        return java.util.Arrays.copyOf(tempArr,n);
    }

    static double[] arraySample(int n, double[] array)
    {
        double[] tempArr = array.clone();
        arrayShuffle(tempArr);
        return java.util.Arrays.copyOf(tempArr,n);
    }

    static int[] arraySampleProb(int n, int[] array, double[] probs)
    {
        int pos;
        double rand;
        int[] newArr = new int[n];
        double[] cumProbs = probs.clone();
        Auxils.arrayCumSum(cumProbs);
        Auxils.arrayDiv(cumProbs, cumProbs[cumProbs.length-1]);
        for (int i = 0; i < n; i++) {
            rand = random.nextDouble();
            pos = arraySearch(cumProbs, rand);
            newArr[i] = array[pos];
        }
        return newArr;
    }

    static double[] arraySampleProb(int n, double[] array, double[] probs)
    {
        int pos;
        double rand;
        double[] newArr = new double[n];
        double[] cumProbs = probs.clone();
        Auxils.arrayCumSum(cumProbs);
        Auxils.arrayDiv(cumProbs, cumProbs[cumProbs.length-1]);
        for (int i = 0; i < n; i++) {
            rand = random.nextDouble();
            pos = arraySearch(cumProbs, rand);
            newArr[i] = array[pos];
        }
        return newArr;
    }

    static int randIntProb(int end, double[] probs)
    {
        int val;
        double rand;
        double[] cumProbs = java.util.Arrays.copyOf(probs, end);
        Auxils.arrayCumSum(cumProbs);
        Auxils.arrayDiv(cumProbs, cumProbs[cumProbs.length-1]);
        rand = random.nextDouble();
        val = arraySearch(cumProbs, rand);
        return val;
    }

    static int arraySearch(double[] array, double key)
    {
        int lo = 0;
        int hi = array.length - 1;
        int mid;
        if (key <= array[lo])
            return lo;
        else {
            while (lo < hi) {
                mid = lo + (hi - lo)/2;
                if (key <= array[mid])
                    hi = mid - 1;
                else if (key > array[mid])
                    lo = mid + 1;
            }
            if (key <= array[lo])
                return lo;
            else
                return ++lo;
        }
    }
    
    static int[] enumArray(int from, int to)
    {
        int[] newArr = new int[to-from+1];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = from++;
        return newArr;
    }
    
    static int[] arrayElements(int[] array, int[] pos)
    {
        int[] newArr = new int[pos.length];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = array[pos[i]];
        return newArr;
    }
    
    static boolean[] arrayElements(boolean[] array, int[] pos)
    {
        boolean[] newArr = new boolean[pos.length];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = array[pos[i]];
        return newArr;
    }
    
    static double[] arrayElements(double[] array, int[] pos)
    {
        double[] newArr = new double[pos.length];
        for (int i = 0; i < newArr.length; i++)
            newArr[i] = array[pos[i]];
        return newArr;
    }
    
    static double arrayMean(int[] array)
    {
        double mean = 0;
        for (int i = 0; i < array.length; i++)
            mean += array[i];
        mean /= array.length;
        return mean;
    }

    static double arrayMean(boolean[] array)
    {
        double mean = 0;
        for (int i = 0; i < array.length; i++)
            if (array[i])
                mean++;
        mean /= array.length;
        return mean;
    }

    static double arrayMean(double[] array)
    {
        double mean = 0;
        for (int i = 0; i < array.length; i++)
            mean += array[i];
        mean /= array.length;
        return mean;
    }

    static double arrayMean(int[] array, int end)
    {
        double mean = 0;
        for (int i = 0; i < end; i++)
            mean += array[i];
        mean /= end;
        return mean;
    }

    static double arrayMean(boolean[] array, int end)
    {
        double mean = 0;
        for (int i = 0; i < end; i++)
            if (array[i])
                mean++;
        mean /= end;
        return mean;
    }

    static double arrayMean(double[] array, int end)
    {
        double mean = 0;
        for (int i = 0; i < end; i++)
            mean += array[i];
        mean /= end;
        return mean;
    }

    static int arraySum(int[] array) 
    {
        int sum = 0;
        for (int i = 0; i < array.length; i++)
            sum += array[i];
        return sum;
    }

    static int arraySum(boolean[] array) 
    {
        int sum = 0;
        for (int i = 0; i < array.length; i++)
            if (array[i])
                sum++;
        return sum;
    }

    static double arraySum(double[] array) 
    {
        double sum = 0;
        for (int i = 0; i < array.length; i++)
            sum += array[i];
        return sum;
    }

    static int arraySum(int[] array, int end) 
    {
        int sum = 0;
        for (int i = 0; i < end; i++)
            sum += array[i];
        return sum;
    }

    static int arraySum(boolean[] array, int end) 
    {
        int sum = 0;
        for (int i = 0; i < end; i++)
            if (array[i])
                sum++;
        return sum;
    }

    static double arraySum(double[] array, int end) 
    {
        double sum = 0;
        for (int i = 0; i < end; i++)
            sum += array[i];
        return sum;
    }

    static void arrayCumSum(int[] array) 
    {
            for (int i = 1; i < array.length; i++)
                array[i] += array[i-1];
    }

    static void arrayCumSum(double[] array) 
    {
            for (int i = 1; i < array.length; i++)
                array[i] += array[i-1];
    }

    static void arrayAdd(int[] array, int a) 
    {
        for (int i = 0; i < array.length; i++)
            array[i] += a;
    }

    static void arrayAdd(double[] array, double a) 
    {
        for (int i = 0; i < array.length; i++)
            array[i] += a;
    }

    static void arrayMult(int[] array, int a) 
    {
        for (int i = 0; i < array.length; i++)
            array[i] *= a;
    }

    static void arrayMult(double[] array, double a) 
    {
        for (int i = 0; i < array.length; i++)
            array[i] *= a;
    }

    static void arrayDiv(int[] array, int a) 
    {
        for (int i = 0; i < array.length; i++)
            array[i] /= a;
    }

    static void arrayDiv(double[] array, double a) 
    {
        for (int i = 0; i < array.length; i++)
            array[i] /= a;
    }

    static int[] arrayConcat(int[] first, int[] second) 
    {
        int[] result = java.util.Arrays.copyOf(first, first.length + second.length);
        System.arraycopy(second, 0, result, first.length, second.length);
        return result;
    }    
    
    static double[] arrayConcat(double[] first, double[] second) 
    {
        double[] result = java.util.Arrays.copyOf(first, first.length + second.length);
        System.arraycopy(second, 0, result, first.length, second.length);
        return result;
    }    

}


