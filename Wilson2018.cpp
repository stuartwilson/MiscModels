/*
 Implementation of a Boolean network model of the evolution of cortical development
 as described by Wilson, Whiteley, and Krubitzer in a manuscript submitted to
 PLoS Computational Biology in 2018.
 */


#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include <ctime>

using namespace std;

vector<vector<int> > graph(vector<vector<int> >);
double randDouble(void);
bool randBool(double);
int bin2int(vector<bool>);
vector<bool> int2bin(long int, long int);


/*
 DEFINE BOOLEAN NETWORK
 */

class Net{

public:
    long int nG, nS, nT;
    vector<bool> G;                 // genome
    vector<vector<int> > A, C;      // attractors, connections
    vector<int> I;                  // attractor ID
    double normVal;
    long int initA, initB, targetA, targetB;

    Net(long int nG){

        this->nG = nG;          // number of genes
        nS = pow(2,nG);         // number of states
        nT = pow(2,nG-1);       // number of truth table combinations
        normVal = 1./(pow(2,nG+1)); // for fitness function

        C.resize(nS);
        for(int s=0;s<nS;s++){ C[s].resize(2); }

        G.resize(nG*nT);
        I.resize(nS);

        // define initial and target states
        initA = 0;
        initB = pow(2,nG-1);
        vector<bool> tarA(nG,false);
        vector<bool> tarB(nG,false);
        for(int i=0;i<nG;i++){
            tarA[i]=i%2;
            tarB[i]=(i+1)%2;
        }
        targetA = bin2int(tarA);
        targetB = bin2int(tarB);

        randomGenome();
    }

    // generate a random genome
    void randomGenome(void){
        for(long int i=0;i<nG*nT;i++){
            G[i]=randBool(0.5);
        }
    }

    vector<bool> getNext(vector<bool> S){
        vector<bool> Snext(nG,false);
        for(long int i=0;i<nG;i++){
            int k = 0;
            int s = 1;
            for(int j=nG;j--;){
                if(i!=j){
                    k+=s*S[j];
                    s*=2;
                }
            }
            Snext[i] = G[i*nT+k];
        }
        return Snext;
    }

    // return number of iterations before a state is repeated
    int getDist(long int s){
        int d = 0;
        vector<bool> V (nS,false); // visited states
        vector<bool> S = int2bin(s,nG);
        while(true){
            int v = bin2int(S);
            if(V[v]) break;
            V[v] = true;
            d++;
            S = getNext(S);
        }
        return d-1;
    }


    // return whether there is a path from a to b
    bool isPath(long int start, long int finish){
        vector<bool> S = int2bin(start,nG);
        for(int i=0;i<=nS;i++){
            S = getNext(S);
            if (bin2int(S)==finish) return true;
        } return false;
    }


    // compute the attractor landscape
    void getAttractors(void){

        for(int s=0;s<nS;s++){
            C[s][0] = s;
            C[s][1] = bin2int(getNext(int2bin(s,nG)));
        }
        A = graph(C);

        for(int i=0;i<A.size();i++){
            for(int j=0;j<A[i].size();j++){
                I[A[i][j]] = i;
            }
        }
    }

    double fitness1(void){ // force to be on separate attractors
        getAttractors();
        if( (I[initA]==I[targetA])
            && (I[initB]==I[targetB])
            && (I[targetA]!=I[targetB])
            ){

            int DtaA = getDist(targetA);
            int DtaB = getDist(targetB);
            if ((getDist(initA)>=DtaA) && (getDist(initB)>=DtaB)){
                return 1.0-(DtaA+DtaB)*normVal;
            }
        }
        return 0.;
    }

    double fitness2(void){ // don't force on separate attractors
        getAttractors();
        if( (I[initA]==I[targetA])
            && (I[initB]==I[targetB])
            ){

            int DtaA = getDist(targetA);
            int DtaB = getDist(targetB);
            if ((getDist(initA)>=DtaA) && (getDist(initB)>=DtaB)){
                return 1.0-(DtaA+DtaB)*normVal;
            }
        }
        return 0.;
    }

    double fitness3(void){ // don't force on separate attractors

        if (isPath(initA,targetA) && isPath(initB,targetB)){
            int DtaA = getDist(targetA);
            int DtaB = getDist(targetB);
            if ((getDist(initA)>=DtaA) && (getDist(initB)>=DtaB)){
                return 1.0-(DtaA+DtaB)*normVal;
            } else {
                return 0.;
            }
        } else {
            return 0.;
        }


    }


    double fitness4(void){ // don't force on separate attractors

        if (isPath(initA,targetA) &&
            isPath(initB,targetB) &&
            isPath(targetA,targetA) &&
            isPath(targetB,targetB)){
            int DtaA = getDist(targetA);
            int DtaB = getDist(targetB);
            if ((getDist(initA)>=DtaA) && (getDist(initB)>=DtaB)){
                return 1.0-(DtaA+DtaB)*normVal;
            } else {
                return 0.;
            }
        } else {
            return 0.;
        }


    }

};


// MAIN EVOLUTIONARY LOOP

int main(int argc, char **argv){

    if(argc<6){
        cout<<"useage: ./model filename genes mutation seed generations"<<endl<<"e.g.  : ./model output 5 0.05 1"<<endl;
        return 0;
    }
    // ./model file Ngenes Prob Seed nGens

    ofstream outFile; std::stringstream oFile; oFile<<argv[1]<<".bin";
    outFile.open(oFile.str().c_str(),ios::out|ios::binary);

    int nGenes = atoi(argv[2]);
    double mutationProbability = atof(argv[3]);
    int seed = atoi(argv[4]);
    int nGenerations = atoi(argv[5]);

    Net N(nGenes);

    int interval = 0;
    double fa = 0.;

    int gens25 = nGenerations/4;
    int gens50 = gens25*2;
    int gens75 = gens25*3;

    for (int t=0;t<nGenerations;t++){

        if(t==gens25) cout<<endl<<" 25% "<<endl;
        if(t==gens50) cout<<endl<<" 50% "<<endl;
        if(t==gens75) cout<<endl<<" 75% "<<endl;

        if (fa == 1.0){
            cout<<"."<<flush;
            outFile.write((char*)(&interval),sizeof(int));
            interval = 0;
            N.randomGenome();
            fa = 0.;

        } else {
            vector<bool> Gpre = N.G;
            for(int i=0;i<Gpre.size();i++){
                if(randBool(mutationProbability)) N.G[i]=~N.G[i];
            }
            double fb = N.fitness4();

            if(fb<fa) N.G=Gpre;
            else fa=fb;

            interval++;
        }
    }
    outFile.close();
    return 0;
}







/*
 FUNCTIONS
 */

// RETURN A RANDOM DOUBLE BETWEEN 0 AND 1
double randDouble(void){
    return ((double) rand())/(double)RAND_MAX;
}

// RETURN A BOOL WITH PROBABILITY OF BEING TRUE GIVEN BY P
bool randBool(double P){
    return (randDouble()<P);
}

// CONVERT A BINARY (VECTOR OF BOOLS) TO AN INTEGER
int bin2int(vector<bool> b){
    int k = 0;
    int s = 1;
    for(int i=b.size();i--;){
        k += s*b[i];
        s *= 2;
    }
    return k;
}

// CONVERT AN INTEGER TO A BINARY (VECTOR OF BOOLS)
vector<bool> int2bin(long int i, long int n){
    vector<bool> x(n,false);
    long int a = pow(2,n-1);
    for(long int k=0;k<n;k++){
        x[k] = (i>=a);
        if(x[k]) i-=a;
        a/=2;
    }
    return x;
}


// TAKES A LIST OF PAIRS AND RETURNS A VECTOR OF CLUSTERS SORTED BY CLUSTER SIZE
vector<vector<int> > graph(vector<vector<int> > conn){
    int N = conn.size();
    vector<vector<int> > graph;
    int a, b, Ain, Bin;
    bool nov, ain, bin;
    for(int k=0;k<conn.size();k++){
        a = conn[k][0];
        b = conn[k][1];
        if(a != b){
            // reset flags
            Ain = -1;
            Bin = -1;
            nov = true;
            // check if conn k is in the graph
            int I = graph.size();
            for(int i=0;i<I;i++){
                ain = false;
                bin = false;
                int J = graph[i].size();
                for(int j=0;j<J;j++){
                    if(a == graph[i][j]){
                        ain = true;
                        Ain = i;
                    }
                    if(b == graph[i][j]){
                        bin = true;
                        Bin = i;
                    }
                }
                if(ain && !bin) graph[i].push_back(b);
                if(!ain && bin) graph[i].push_back(a);
                if(ain||bin) nov=false;
            }
            // Add a new group
            if(nov) graph.push_back(conn[k]);

            // Join two existing groups
            if(Ain>-1 && Bin>-1 && Ain!=Bin){
                graph[Ain].pop_back();
                graph[Bin].pop_back();
                for(int l=0;l<graph[Bin].size();l++){
                    graph[Ain].push_back(graph[Bin][l]);
                }
                graph.erase(graph.begin()+Bin);
            }
        }
    }

    for(int k=0;k<N;k++){
        bool isolated = true;
        int I = graph.size();
        for(int i=0;i<I;i++){
            int J = graph[i].size();
            for(int j=0;j<J;j++){
                if(k==graph[i][j]){
                    isolated = false;
                    //break;
                }
            }
        }
        if(isolated){
            vector<int> isolate(1,k);
            graph.push_back(isolate);
        }
    }

    return graph;
}



