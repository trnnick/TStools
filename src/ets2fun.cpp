#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
//using namespace arma;
 
arma::mat matrixpower(arma::mat A, int power){
    arma::mat B = A;
    if(power>1){
        for(int i = 1; i < power; i=i+1){
            B = B * A;
        }
    }
    else if(power==0){
        B.eye();
    }
    return B;
}

double errorf(double yact, double yfit, char Etype){
    if(Etype=='A'){
        return yact - yfit;
    }
    else{
        return (yact - yfit)/yfit;
    }
}

/* Function will be needed to estimate the correct error for ETS when trace model selection with r(xt) is sorted out. */
arma::mat errorvf(arma::mat yact, arma::mat yfit, char Etype){
    if(Etype=='A'){
        return yact - yfit;
    }
    else{
        return (yact - yfit) / yfit;
    }
}

arma::rowvec rvalue(arma::rowvec xt, arma::rowvec matw, char Etype, char Ttype, char Stype, int ncomponents){
    arma::rowvec r(ncomponents);
    arma::rowvec xtnew;
    arma::vec matwnew = arma::vec(matw.memptr(),ncomponents,false,true);
    
    xtnew = xt % matw;

    if(Etype=='A'){
        if(Ttype=='N'){
            if(Stype!='M'){
                r.ones();
            }
            else{
                r(0) = 1 / sum(xtnew(arma::span(1, ncomponents-1)));
                r(arma::span(1, ncomponents-1)).fill(1 / xtnew(0));
            }
        }
        else if(Ttype=='A'){
            if(Stype!='M'){
                r.ones();
            }
            else{
                r(0) = 1 / sum(xtnew(arma::span(2, ncomponents-1)));
                r(1) = 1 / sum(xtnew(arma::span(2, ncomponents-1)));
                r(arma::span(2, ncomponents-1)).fill(1 / sum(xtnew(arma::span(0,1))));
            }
        }
        else if(Ttype=='M'){
            if(Stype=='N'){
                r(0) = 1;
                r(1) = 1 / xtnew(0);
            }
            else if(Stype=='A'){
                r.ones();
                r(1) = 1 / xtnew(0);
            }
            else if(Stype=='M'){
                r(0) = 1 / sum(xtnew(arma::span(2, ncomponents-1)));
                r(1) = 1 / (xtnew(0) * sum(xtnew(arma::span(2, ncomponents-1))));
                r(arma::span(2, ncomponents-1)).fill(1 / (xt(0) * pow(xt(1),matw(2))));
            }
        }
    }
    else{
        if(Ttype=='N'){
            if(Stype=='N'){
                r(0) = xtnew(0);
            }
            else if(Stype=='A'){
                r.fill(sum(xt * matwnew));
            }
            else if(Stype=='M'){
                r(0) = xt(0);
                r(arma::span(1, ncomponents-1)) = xt(arma::span(1, ncomponents-1));
            }
        }
        else if(Ttype=='A'){
            if(Stype!='M'){
                r.fill(sum(xt * matwnew));
            }
            else if(Stype=='M'){
                r(0) = sum(xt(arma::span(0,1)) * matwnew(arma::span(0,1)));
                r(1) = r(0);
                r(arma::span(2, ncomponents-1)) = xt(arma::span(2, ncomponents-1));
            }
        }
        else if(Ttype=='M'){
            if(Stype=='N'){
                r(0) = exp(sum(log(xt(arma::span(0,1))) * matwnew(arma::span(0,1))));
                r(1) = exp(matw(1) * log(xt(1)));
            }
            else if(Stype=='A'){
                r(0) = exp(sum(log(xt(arma::span(0,1))) * matwnew(arma::span(0,1)))) + sum(xtnew(arma::span(2, ncomponents-1)));
                r(1) = r(0) / xt(0);
                r(arma::span(2, ncomponents-1)).fill(r(0));
            }
            else if(Stype=='M'){
                r(0) = exp(sum(log(xt(arma::span(0,1))) * matwnew(arma::span(0,1))));
                r(1) = exp(matw(1) * log(xt(1)));
                r(arma::span(2, ncomponents-1)).fill(sum(xtnew(arma::span(2, ncomponents-1))));
            }
        }
    }

    return r;
}

// [[Rcpp::export]]
RcppExport SEXP fitets2(SEXP xt, SEXP F, SEXP w, SEXP yt, SEXP g, SEXP Etype, SEXP Ttype, SEXP Stype, SEXP sf) {
    NumericMatrix mxt(xt);
    NumericMatrix mF(F);
    NumericMatrix vw(w);
    NumericMatrix vyt(yt);
    NumericMatrix vg(g);
    arma::mat matxt(mxt.begin(), mxt.nrow(), mxt.ncol(), false);
    arma::mat matF(mF.begin(), mF.nrow(), mF.ncol(), false);
    arma::mat matw(vw.begin(), vw.nrow(), vw.ncol(), false);
    arma::mat matyt(vyt.begin(), vyt.nrow(), vyt.ncol(), false);
    arma::mat matg(vg.begin(), vg.nrow(), vg.ncol(), false);
    int freq = as<int>(sf);
    char E = as<char>(Etype);
    char T = as<char>(Ttype);
    char S = as<char>(Stype);
    int obs = matyt.n_rows;
    int freqtail;
    int j;
    int ncomponents;
    int ncomponentsall;

    arma::mat matyfit(obs, 1, arma::fill::zeros);
    arma::mat materrors(obs, 1, arma::fill::zeros);

    arma::mat matxtnew;
    arma::mat matwnew;
    arma::mat matFnew;
    arma::mat matgnew;
    arma::mat dummy(freq,freq, arma::fill::eye);

/* xt is transformed into obs by n.components+freq matrix, where last freq are seasonal coefficients,
w is matrix obs by n.components+freq with dummies in freq parts,
F is a matrix as in Hyndman et al.
g is the vector of 
dummy contains dummy variables for seasonal coefficients
*/
    if(S!='N'){
        ncomponents = mxt.ncol()-1;
        ncomponentsall = ncomponents + freq;

        matgnew.set_size(obs,ncomponentsall);
        matgnew.cols(0,ncomponents-1).each_row() = trans(matg.rows(0,ncomponents-1));

/* Fill in xt with provided initial level and trend and then - with the provided initial seasonals */
        matxtnew.set_size(obs+1,ncomponentsall);
        matxtnew.zeros();
        matxtnew.submat(0, 0, 0, ncomponents-1) = matxt.submat(0, 0, 0, ncomponents-1);
        matxtnew.submat(0, ncomponents, 0, ncomponentsall-1) = trans(matxt.submat(0, ncomponents, freq-1, ncomponents));

        matwnew.set_size(obs,ncomponentsall);
        matFnew.eye(ncomponentsall,ncomponentsall);
/* Fill in the matrix w and cube F with the provided values of w and F */
        for(int i=0; i<ncomponents; i=i+1){
            matwnew.col(i).each_row() = matw.col(i);
            matFnew(i,0) = matF(i,0);
            if(ncomponents>1){
                matFnew(i,1) = matF(i,1);
            }
        }

/* Fill in dummies for the seasonal components */
        for(int i=0; i < std::floor(obs/freq); i=i+1){
            matwnew.submat(i*freq, ncomponents, (i+1)*freq-1, ncomponentsall-1) = dummy;
            matgnew.submat(i*freq, ncomponents, (i+1)*freq-1, ncomponentsall-1) = dummy * matg(ncomponents,0);
        }

        freqtail = obs - freq * std::floor(obs/freq);
        if(freqtail!=0){
            j = 0;
            for(int i=obs-freqtail; i<obs; i=i+1){
                matwnew.submat(i, ncomponents, i, ncomponents-1+freq) = dummy.row(j);
                matgnew.submat(i, ncomponents, i, ncomponents-1+freq) = dummy.row(j) * matg(ncomponents,0);
                j=j+1;
            }
        }
    }
    else{
        ncomponents = mxt.ncol();
        ncomponentsall = ncomponents;
        matxtnew = matxt;
        matgnew.set_size(obs,ncomponents);
        matgnew.each_row() = trans(matg);

        matwnew.set_size(obs,ncomponents);
        matFnew.eye(ncomponents,ncomponents);

/* Fill in the matrix w and cube F with the provided values of w and F */
        for(int i=0; i<ncomponents; i=i+1){
            matwnew.col(i).each_row() = matw.col(i);
            matFnew(i,0) = matF(i,0);
            if(ncomponents>1){
                matFnew(i,1) = matF(i,1);
            }
        }
    }

    if(T!='M' & S!='M'){
        for (int i=0; i<obs; i=i+1) {
            matyfit.row(i) = matwnew.row(i) * trans(matxtnew.row(i));
            materrors(i,0) = errorf(matyt(i,0), matyfit(i,0), E);
            matxtnew.row(i+1) = matxtnew.row(i) * trans(matFnew) + trans(materrors.row(i)) * matgnew.row(i) % rvalue(matxtnew.row(i), matwnew.row(i), E, T, S, ncomponentsall);
        }
    }
    else if(T!='A' & S!='A'){
        if(T=='M' | S=='M'){
            for (int i=0; i<obs; i=i+1) {
                matyfit.row(i) = exp(matwnew.row(i) * log(trans(matxtnew.row(i))));
                materrors(i,0) = errorf(matyt(i,0), matyfit(i,0), E);
                matxtnew.row(i+1) = exp(log(matxtnew.row(i)) * trans(matFnew)) + trans(materrors.row(i)) * matgnew.row(i) % rvalue(matxtnew.row(i), matwnew.row(i), E, T, S, ncomponentsall);
                matxtnew.elem(find(matxtnew.row(i+1)<0)).zeros();
            }
        }
    }
    else if(T=='A' & S=='M'){
        for (int i=0; i<obs; i=i+1) {
            matyfit.row(i) = matwnew.row(i).cols(0,1) * trans(matxtnew.row(i).cols(0,1)) * sum(matwnew.row(i).cols(2,ncomponentsall-1) % matxtnew.row(i).cols(2,ncomponentsall-1));
            materrors(i,0) = errorf(matyt(i,0), matyfit(i,0), E);
            matxtnew.row(i+1) = matxtnew.row(i) * trans(matFnew) + trans(materrors.row(i)) * matgnew.row(i) % rvalue(matxtnew.row(i), matwnew.row(i), E, T, S, ncomponentsall);
            matxtnew.elem(find(matxtnew.row(i+1).cols(2,ncomponentsall-1)<0)).zeros();
        }
    }
    else if(T=='M' & S=='A'){
        for (int i=0; i<obs; i=i+1) {
            matyfit.row(i) = exp(matwnew.row(i).cols(0,1) * log(trans(matxtnew.row(i).cols(0,1)))) + sum(matwnew.row(i).cols(2,ncomponentsall-1) % matxtnew.row(i).cols(2,ncomponentsall-1));
            materrors(i,0) = errorf(matyt(i,0), matyfit(i,0), E);
            matxtnew.row(i+1).cols(0,1) = exp(log(matxtnew.row(i).cols(0,1)) * trans(matFnew.submat(0,0,1,1))) + trans(materrors.row(i)) * matgnew.row(i).cols(0,1) % rvalue(matxtnew.row(i), matwnew.row(i), E, T, S, ncomponentsall).cols(0,1);
            if(arma::sum(matxtnew(i+1,1))<0){
                matxtnew(i+1,1) = matxtnew(i,1);
            }
            matxtnew.row(i+1).cols(2,ncomponentsall-1) = matxtnew.row(i).cols(2,ncomponentsall-1) * trans(matFnew.submat(2,2,ncomponentsall-1,ncomponentsall-1)) + trans(materrors.row(i)) * matgnew.row(i).cols(2,ncomponentsall-1) % rvalue(matxtnew.row(i), matwnew.row(i), E, T, S, ncomponentsall).cols(2,ncomponentsall-1);
        }
    }

// Fill in xt for the seasonal models
    if(S!='N'){
        matxt.submat(freq-1,0,mxt.nrow()-1,ncomponents-1) = matxtnew.cols(0,ncomponents-1);

        for(int i=0; i < std::floor(obs/freq); i=i+1){
            matxt.submat((i+1)*freq,ncomponents,(i+2)*freq-1,ncomponents) = matxtnew.submat(i*freq+1,ncomponents,(i+1)*freq,ncomponentsall-1).diag();
        }

        if(freqtail!=0){
            j = 0;
            for(int i=matxtnew.n_rows-freqtail; i<matxtnew.n_rows; i=i+1){
                matxt(i+freq-1, ncomponents) = matxtnew(i, ncomponents+j);
                j=j+1;
            }
        }
    }
    else{
        matxt = matxtnew;
    }

    return List::create(Named("xt") = matxt, Named("yfit") = matyfit, Named("errors") = materrors);
}

// [[Rcpp::export]]
RcppExport arma::mat forets2(SEXP xt, SEXP F, SEXP w, SEXP h, SEXP Ttype, SEXP Stype, SEXP sf) {
    NumericMatrix mxt(xt);
    NumericMatrix mF(F);
    NumericMatrix vw(w);
    arma::mat matxt(mxt.begin(), mxt.nrow(), mxt.ncol(), false);
    arma::mat matF(mF.begin(), mF.nrow(), mF.ncol(), false);
    arma::mat matw(vw.begin(), vw.nrow(), vw.ncol(), false);
    int hor = as<int>(h);
    int freq = as<int>(sf);
    char T = as<char>(Ttype);
    char S = as<char>(Stype);
    int hh;

    arma::mat matyfor(hor, 1, arma::fill::zeros);
    arma::mat matxtnew;
    arma::mat matwnew;
    arma::mat matFnew;
    arma::mat seasxt(hor, 1, arma::fill::zeros);

    if(S!='N'){
        hh = std::min(hor,freq);
        seasxt.submat(0, 0, hh-1, 0) = matxt.submat(0,matxt.n_cols-1,hh-1,matxt.n_cols-1);

        if(hor > freq){
            for(int i = freq; i < hor; i=i+1){
                seasxt.row(i) = seasxt.row(i-freq);
            }
        }

        matxtnew = matxt.submat(mxt.nrow()-1, 0, mxt.nrow()-1, mxt.ncol()-2);
        matwnew = matw.cols(0, vw.ncol()-2);
        matFnew = matF.submat(0, 0, mF.nrow()-2, mF.ncol()-2);
    }
    else{
        matxtnew = matxt.submat(mxt.nrow()-1, 0, mxt.nrow()-1, mxt.ncol()-1);
        matwnew = matw;
        matFnew = matF;
    }

    if(T!='M'){
        for(int i = 0; i < hor; i=i+1){
            matyfor.row(i) = matwnew * matrixpower(matFnew,i) * trans(matxtnew);
        }
    }
    else{
        for(int i = 0; i < hor; i=i+1){
            matyfor.row(i) = exp(matwnew * matrixpower(matFnew,i) * log(trans(matxtnew)));
        }
    }
    if(S=='A'){
        matyfor = matyfor + seasxt;
    }
    else if(S=='M'){
        matyfor = matyfor % seasxt;
    }

    return matyfor;
}

// [[Rcpp::export]]
RcppExport arma::mat errorets2(SEXP xt, SEXP F, SEXP w, SEXP yt, SEXP h, SEXP Etype, SEXP Ttype, SEXP Stype, SEXP sf, SEXP trace) {
    NumericMatrix mxt(xt);
    NumericMatrix mF(F);
    NumericMatrix vw(w);
    NumericMatrix vyt(yt);
    arma::mat matxt(mxt.begin(), mxt.nrow(), mxt.ncol(), false);
    arma::mat matF(mF.begin(), mF.nrow(), mF.ncol(), false);
    arma::mat matw(vw.begin(), vw.nrow(), vw.ncol(), false);
    arma::mat matyt(vyt.begin(), vyt.nrow(), vyt.ncol(), false);
    int hor = as<int>(h);
    int freq = as<int>(sf);
    char E = as<char>(Etype);
    char T = as<char>(Ttype);
    char S = as<char>(Stype);
    bool tr = as<bool>(trace);
    int obs = vyt.nrow();
    int hh;
    arma::mat materrors;

    if(tr){
        materrors.set_size(obs, hor);
        materrors.fill(NA_REAL);
    }
    else{
        materrors.set_size(obs, 1);
    }

    if(tr==true){
        for(int i = 0; i < obs; i=i+1){
            hh = std::min(hor, obs-i);
/*            materrors.submat(i, 0, i, hh-1) = trans(errorvf(matyt.rows(i, i+hh-1),forets2(wrap(matxt.rows(i,i+freq-1)),F,w,wrap(hh),Ttype,Stype,sf),E)); */
            materrors.submat(i, 0, i, hh-1) = trans(matyt.rows(i, i+hh-1) - forets2(wrap(matxt.rows(i,i+freq-1)),F,w,wrap(hh),Ttype,Stype,sf));
        }
    }
    else{
	    for(int i = 0; i < obs; i=i+1){
/*	        materrors.row(i) = errorvf(matyt.row(i),forets2(wrap(matxt.rows(i,i+freq-1)),F,w,wrap(1),Ttype,Stype,sf),E); */
            materrors.row(i) = matyt.row(i) - forets2(wrap(matxt.rows(i,i+freq-1)),F,w,wrap(1),Ttype,Stype,sf);
	    }
    }

    return materrors;
}