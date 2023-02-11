#include <ilcplex/ilocplex.h>
#include <vector>
#include <fstream>

using namespace std;

#define RC_EPS 1.0e-6

// can be printed in an output.txt file
void output(IloCplex& cutSolver, IloNumVarArray &Cut, vector<vector<int> > &patterns, int iter)
{
    ofstream out("/Users/asc/cplex/cplex/output.txt");
    out << "total stocks required: " << cutSolver.getObjValue() << endl;
    out << endl;
    for (IloInt j = 0; j < Cut.getSize(); j++) {
        
        for(int x=0;x<patterns[0].size();x++){
            out << patterns[j][x] << " ";
        }
        out << "cut" << j << " = " << cutSolver.getValue(Cut[j]) << endl;
    }
    
}

int rollwidth;
vector<double> sz;
vector<int> amt;

void input(){
    ifstream in("/Users/asc/cplex/cplex/input.txt");
    string s;
    getline(in,s);
    rollwidth = stoi(s);
    
    string st;
    getline(in,st);
    st=st+" ";

    string x="";
    for(int i=0;i<st.size();i++){
        if(st[i] == ' '){
            double tmp = stod(x);
            sz.push_back(tmp);
            x="";
        }
        else
            x.push_back(st[i]);
    }

    string str;
    getline(in,str);
    str=str+" ";
    x="";
    for(int i=0;i<str.size();i++){
        if(str[i] == ' '){
            int tmp = stoi(x);
            amt.push_back(tmp);
            x="";
        }
        else
            x.push_back(str[i]);
    }
    return;
}


int main(){
    
    IloEnv env;
    
    //main problem
    IloModel cutOpt (env);
    
    IloObjective StocksUsed = IloMinimize(env);
    cutOpt.add(StocksUsed);
    
    input(); // input to be taken from input.txt file
    
    IloNum rollWidth = rollwidth;
    
    IloNumArray size(env);
    for(int i=0;i<sz.size();i++){
        size.add(sz[i]);
    }
    IloNumArray amount(env);
    for(int i=0;i<amt.size();i++){
        amount.add(amt[i]);
    }
    
    vector< vector<int> > patterns;
    for(int x=0;x<amount.getSize();x++){
        vector<int> tmp(amount.getSize(),0);
        tmp[x]=rollwidth / size[x];
        patterns.push_back(tmp);
    }
    
    // log file
    ofstream log("/Users/asc/cplex/cplex/log.txt");
    
    log << "Initial Patterns added to the master problem" << endl;
    for(auto x: patterns){
        for(auto y: x){
            log<<y<<" ";
        }
        log<<endl;
    }
    log << endl;
//    log << rollWidth <<endl;
//    log << size << endl;
//    log << amount << endl;
    
    IloRangeArray  Fill = IloAdd(cutOpt, IloRangeArray(env, amount, IloInfinity));
    
    IloNumVarArray Cut(env);
    
    IloInt nWdth = size.getSize();
    for (IloInt j = 0; j < nWdth; j++) {
        Cut.add(IloNumVar(StocksUsed(1) + Fill[j](int(rollWidth / size[j]))));
    }
    
    IloCplex cutSolver(cutOpt);
    
    // subproblem
    IloModel patGen (env);

    IloObjective ReducedCost = IloMinimize(env, 1);
    patGen.add(ReducedCost);

    IloNumVarArray Use(env, nWdth, 0.0, IloInfinity, ILOINT);
    patGen.add(IloScalProd(size, Use) <= rollWidth);

    IloCplex patSolver(patGen);
    
    IloNumArray price(env, nWdth);
    IloNumArray newPatt(env, nWdth);
    
    IloNumArray dual(env, nWdth);
    
    int iter=0;
    
    while(1){
        
        log << "Iteration " << ++iter << ":" << endl;
        
        cutSolver.solve();
        
        log << "Objective value of the master-problem: " << cutSolver.getObjValue() << endl;
       
        for (IloInt i = 0;i < nWdth;i++){
            price[i] = -cutSolver.getDual(Fill[i]);
            dual[i] = -price[i];
        }
        log << "Dual values: " << dual << endl;
        
        ReducedCost.setLinearCoefs(Use, price);
        
        patSolver.solve();
        log << "Objective value of the sub-problem: " << patSolver.getObjValue() << endl;
        
        if (patSolver.getValue(ReducedCost) > -RC_EPS) break;
       
        patSolver.getValues(newPatt, Use);
        Cut.add(IloNumVar(StocksUsed(1) + Fill(newPatt)));
        
        vector<int> tmp;
        for(int x=0;x<newPatt.getSize();x++){
            tmp.push_back(newPatt[x]);
        }
        patterns.push_back(tmp);
        
        log << "New Pattern generated: " << newPatt << endl << endl;
    
    }
    cutOpt.add(IloConversion(env, Cut, ILOINT));
    
    cutSolver.solve();
    
    output(cutSolver, Cut, patterns, iter);
    
    env.end();
    
    return 0;
}
