%linux update
%nuweb constrainedLinearPredictor
%chmod +x theChMods
%./theChmods
%./rmTheDirs rid of links
%./mkTheDirs
%./makeLns
%rm nw*
%nuweb constrainedLinearPredictor

%ln -s ATconstrainedLinearPredictor @@constrainedLinearPredictor 

%dvips -o silly.ps -p 150 secondOrderCode
%mpage -t silly.ps| lpr
%dvips -o silly.ps -l 149 secondOrderCode
%mpage -t silly.ps| lpr
%remember to cvs co -P to avoid checking out dirs
%chmod +x doLinks.sh(*needed in unix but not cygwin*)
%./doLinks.sh
%nuweb secondOrderCode.w
%java -jar /msu/res5/jarFiles/argoUml/argouml.jar &
% secondOrderCode.zargo
%addpath('/msu/home/m1gsa00/sp_solve')
\documentclass[12pt]{article} \usepackage{latexsym} \usepackage{verbatim}
\usepackage{moreverb}
%\usepackage{notebook}
%\usepackage{alltt}
\usepackage{graphicx} \usepackage{rotating} \usepackage{authordate1-4}
\usepackage{pseudocode}
\usepackage{amsmath}
\usepackage{pifont}


%\includegraphics[width=0.5cm]{/msu/home/m1gsa00/downloads/mathematica2.ps}

\usepackage{makeidx}
\makeindex
\begin{document}

\title{constrainedLinearPredictor: MATLAB Code for Kalman Filter Initialization Experiments\\
$Revision: 1.6 $ 
}
\author{Gary S. Anderson}


\maketitle
\begin{abstract}
  Alejandro Justiniano and Luca Gueirreri have developed code for estimating
DSGE models. This code uses a Kalman Filter to Estimate the Unobservable
Expectations Variables.  This code facilitates the application of
a number of techniques for initializing the Kalman Filter in that context.

\end{abstract}
\tableofcontents

\newpage
\centerline{{\large \bf Things To Do}}

\begin{description}
\item[update distspec] when specifying priors make sure seed other
distribution slots get correct value
\item[optional value] make value optional scalar or interval for initialization
\item[optional transform function] user supplied functions possible
\item[why gamma draw failures?] so many for some settings
how to handle fails in genDraws? reconcile 1e-4 with gamrnd output for
small variance
\item[check how seed works] need to implement seed by distribution?
\end{description}
\newpage

\section{Overview of Calculations}
\label{sec:overview}

\section{setupPath}

@o setupPath.m
@{
function setupPath
addpath ../matWithDrvs;
@}
\section{sKron}

@o sKron.m
@{
function theRes=sKron(aa,bb)
theRes=sparse(kron(aa,bb);
@}
\section{nRows}

@o nRows.m
@{
function theRes=nRows(aMat)
theRes=size(aMat,1);
@}
\section{nCols}

@o nCols.m
@{
function theRes=nCols(aMat)
theRes=size(aMat,2);
@}
\section{isSq}

@o isSq.m
@{
function theRes=isSq(aMat)
theRes=nRows(aMat)==nCols(aMat);
@}
\section{cmpFtt}

@o cmpFtt.m
@{
function theRes=cmpFtt(ff,tt)
theRes=[];
if(isSq(ff));
fDim=nRows(ff);
theRes=zeros(fDim,tt*fDim);
theRes(:,(tt-1)*fDim+(1:fDim))=eye(fDim);
for ii=2:tt
theRes(:,(1:fDim)+(tt-ii)*fDim)=ff*theRes(:,(1:fDim)+(tt-ii+1)*fDim);
end
else
error('cmpFtt: expects square matrix input');
end
@}
\section{cmpGtt}

@o cmpGtt.m
@{
function theRes=cmpGtt(ff,hh,tt)
theRes=[];
if(and(isSq(ff),nCols(ff)==nCols(hh)))
fDim=nRows(ff);
hDim=nRows(hh);
theRes=zeros(tt*hDim,fDim);
theRes((1:hDim),:)=hh * ff;
for ii=2:tt
theRes((1:hDim)+(ii-1)*hDim,:)=theRes((1:hDim)+(ii-2)*hDim,:)*ff;
end
else
error('cmpFtt: expects square ff and conformable hh');
end
@}
\section{cmpBigECov}

@o cmpBigECov.m
@{
function [theVss,theBigU,theWs]=cmpBigECov(ff,vv,hh,ww,tt)
fDim=nRows(ff);
theEnd=fDim*tt;
theVss=zeros(fDim,theEnd);
theWs=zeros(theEnd);
theBigU=zeros(theEnd);
theBigU(1:fDim,(1:fDim))=vv;
theVss(1:fDim,(1:fDim))=vv;
for ii=2:tt
theBigU(1:fDim,(1:fDim)+(ii-1)*fDim)=theBigU(1:fDim,(1:fDim)+(ii-2)*fDim)*ff';
end
theBigU((fDim+1):theEnd,(1:fDim))=...
theBigU((1:fDim),(fDim+1):theEnd)';
theWs((1:fDim),:)=theBigU((1:fDim),:);

theRefMats=theBigU(1:fDim,:);
for ii=2:tt
theRefMats=prePostMult(ff,theRefMats);
theBigU((1:fDim)+(ii-1)*fDim,((ii-1)*fDim+1):theEnd)=...
theBigU((1:fDim)+(ii-2)*fDim,((ii-2)*fDim+1):(theEnd-fDim))+theRefMats;
theVss(:,(1:fDim)+(ii-1)*fDim)=...
theBigU((1:fDim)+(ii-1)*fDim,(1:fDim)+(ii-1)*fDim);
theBigU(((ii)*fDim+1):theEnd,(1:fDim)+(ii-1)*fDim)=...
theBigU((1:fDim)+(ii-1)*fDim,((ii)*fDim+1):theEnd)';
theWs((1:fDim)+(ii-1)*fDim,:)=theBigU((1:fDim)+(ii-1)*fDim,:);
end
hMultMat=kronILeft(hh,tt);
theBigU=hMultMat*theBigU*hMultMat'+kronILeft(ww,tt);
theWs=theWs*hMultMat';
function theMult=prePostMult(ff,theMat)
fDim=nRows(ff);matDim=nCols(theMat);
postMat=kronILeft(ff',(matDim/fDim));
theMult=ff*theMat*postMat(:,1:end-fDim);
@}
\section{cmpECov}

@o cmpECov.m
@{
function theRes=cmpECov(ff,vv,tt,varargin)
theRes=[];
if(and(isSq(vv),and(isSq(ff),nCols(vv)==nCols(ff))))
if(nargin>4)
error('cmpECov: too many args')
end
fDim=size(ff,1);
%ft=cmpFtt(ff,tt);
if(nargin==3)
ss=tt;
else
ss=varargin{1};
end
minTtSs=min(tt,ss);

if(tt>ss)
fLeft=ff^(tt-ss);
fRight=eye(fDim);
elseif(tt<ss)
fLeft=eye(fDim);
fRight=ff^(ss-tt);
else
fLeft=eye(fDim);
fRight=eye(fDim);
end
theRes=zeros(fDim);
for ii=1:minTtSs
theRes=theRes+fLeft*vv*fRight';
fLeft=fLeft*ff;
fRight=fRight*ff;
end
else
error('cmpECov: expects square ff,vv conformable');
end
@}
\section{oldCmpECov}

@o oldCmpECov.m
@{
function theRes=cmpECov(ff,vv,tt,varargin)
theRes=[];
if(and(isSq(vv),and(isSq(ff),nCols(vv)==nCols(ff))))
if(nargin>4)
error('cmpECov: too many args')
end
ft=cmpFtt(ff,tt);
if(nargin==3)
theRes=ft * kron(eye(tt),vv) * ft';
else
ss=varargin{1};
fs=cmpFtt(ff,ss);
theMiddle=cmpMiddle(vv,nRows(ff),tt,ss);
theRes=ft * theMiddle * fs';
end
else
error('cmpECov: expects square ff,vv conformable');
end
if(nargin>3)
ss=varargin{1};
else
ss=tt;
end
theDiff=(theRes-cmpECov(ff,vv,tt,ss));
disp(sprintf('diff in cmpECov=%d(%d,%d)',norm(theDiff),ss,tt));

function theMiddle=cmpMiddle(vv,fRows,tt,ss)
switch sign(tt-ss)
case 1
theMiddle=[kron(eye(ss),vv);zeros((tt-ss)*fRows,ss*fRows)];
case -1
theMiddle=[kron(eye(tt),vv) zeros(tt*fRows,(ss-tt)*fRows)];
case 0
theMiddle=kron(eye(tt),vv);
end

@}
\section{newCmpDelCov}

@o newCmpDelCov.m
@{
function theRes=newCmpDelCov(ff,hh,vv,ww,vtt,tt,varargin)
theRes=[];
if(kfMatsQ(ff,vv,hh,ww))
if(nargin>7)
error('cmpDelCov: wrong number of args')
end
if(nargin==6)
%theRes=hh * vtt * hh' + ww; not used erroneous anyway and about to replace
else
ss=varargin{1};
vts=cmpECov(ff,vv,tt,ss);
theRes=hh * vts *hh';
if(tt==ss)
theRes=theRes+ww;
end
end
else
error('cmpDelCov: problem with kfMatsQ')
end

@}
\section{cmpDelCov}

@o cmpDelCov.m
@{
function theRes=cmpDelCov(ff,hh,vv,ww,tt,varargin)
theRes=[];
if(kfMatsQ(ff,vv,hh,ww))
if(nargin>6)
error('cmpDelCov: wrong number of args')
end
vtt=cmpECov(ff,vv,tt);
if(nargin==5)
theRes=hh * vtt * hh' + ww;
else
ss=varargin{1};
vts=cmpECov(ff,vv,tt,ss);
theRes=hh * vts *hh';
if(tt==ss)
theRes=theRes+ww;
end
end
else
error('cmpDelCov: problem with kfMatsQ')
end

@}
\section{kfMatsQ}

@o kfMatsQ.m
@{
function theRes=kfMatsQ(ff,vv,hh,ww)
theRes=1;
theRes=and(isSq(ww),and(isSq(ff),isSq(vv)));
if(and(~isempty(vv),~isempty(ff)))
theRes=and(theRes,nRows(vv)==nRows(ff));
end
if(and(~isempty(hh),~isempty(ff)))
theRes=and(theRes,nCols(hh)==nRows(ff));
end
if(and(~isempty(hh),~isempty(ww)))
theRes=and(theRes,nRows(hh)==nRows(ww));
end
@}
\section{cmpEDelCov}

@o cmpEDelCov.m
@{

function theRes=cmpEDelCov(ff,hh,vv,tt,ss)
theRes=[];
%if(kfMatsQ(ff,vv,hh,ww))
vts=cmpECov(ff,vv,tt,ss);
theRes=vts * hh';
end
@}
\section{cmpPredErrs}

@o cmpPredErrs.m
@{
function theErrs=cmpPredErrs(ff,hh,vv,ww,tt,yObs)
[xpred,ypred]=cmpXHatsYHats(ff,hh,vv,ww,tt,yObs);
hDim=size(hh,1);
theLen=size(yObs,1);
numObs=theLen/hDim;
theErrs=yObs(1:theLen-hDim)-vec(ypred);
theErrs=reshape(theErrs,hDim,numObs-1);


@}
\section{cmpXHatsYHats}

@o cmpXHatsYHats.m
@{
function [xHats,yHats,theDisps,theErrs,logLH]=cmpXHatsYHats(ff,hh,vv,ww,tt,yObs)
fDim=nRows(ff);
hDim=nRows(hh);
gt=cmpGtt(ff,hh,tt);
%uu=cmpBigU(ff,hh,vv,ww,tt);
[newvss,uu,theWs]=cmpBigECov(ff,vv,hh,ww,tt);
uInv=inv(uu);
fProd=eye(fDim);
xHats=zeros(fDim,tt);
logLH=zeros(1,tt);
theDisps=zeros(hDim,hDim*tt);
yObsRel=yObs(1:(nRows(ww)*tt));
for ss=1:tt
fProd=fProd*ff;
lilU=theWs((1:fDim)+(ss-1)*fDim,:)';%cmpLilU(ff,hh,vv,tt,ss);
bigC=cmpBigC(fProd,gt,lilU,uInv);
xHats(:,ss)=cmpXHatNoInfo(fProd,gt,lilU,uInv,yObsRel,bigC);
vss=newvss(:,(1:fDim)+(ss-1)*fDim);%cmpECov(ff,vv,ss);
[ig,thisDisp]=cmpDispNoInfo(...
fProd,gt,lilU,uInv,vss,bigC);
theDisps(:,(1:hDim)+(ss-1)*hDim)=hh*thisDisp*hh'+ww;
yHats(:,ss)=hh*xHats(:,ss);
theErrs(:,ss)=yObs((1:hDim)+(ss-1)*hDim)-yHats(:,ss);
logLH(1,ss)=onePeriod(theErrs(:,ss),theDisps(:,(1:hDim)+(ss-1)*hDim));
end

%theRes23=cmpDispNoInfo(fs,gt,lilU,uu,vss);


function logLike=onePeriod(theDevs,theVar)
logLike=gaussian_prob(theDevs,zeros(size(theVar,1),1),theVar,true);


@}
\section{priorMeanVar}

@o priorMeanVar.m
@{
function [theRes11,theRes12,theRes21,theRes22,theRes23]=...
priorMeanVar(xZero,pMat,ff,hh,vv,ww,tt,ss,yObs)
fs=ff^ss;
gt=cmpGtt(ff,hh,tt);
yObsRel=yObs(1:(nRows(ww)*tt));
vss=cmpECov(ff,vv,ss);
uu=cmpBigU(ff,hh,vv,ww,tt);
lilU=cmpLilU(ff,hh,vv,tt,ss);
rMat=fs * pMat * fs' + vss;
sMat=fs * pMat * gt' + lilU';
tMat=gt * pMat * gt' + uu;
theRes11=cmpXHat(fs,gt,xZero,sMat,tMat,yObsRel);
theRes12=cmpDisp(rMat,sMat,tMat);
theRes21=cmpXHatNoInfo(fs,gt,lilU,uu,yObsRel);
theRes22=cmpXHatEmpirical(fs,gt,lilU,uu,yObsRel);
theRes23=cmpDispNoInfo(fs,gt,lilU,uu,vss);

@}


\section{cmpXHat}

@o cmpXHat.m
@{
function theRes=cmpXHat(fs,gt,xZero,sMat,tMat,yObs)    
theRes=fs*xZero+sMat* (tMat\ (yObs- gt*xZero));


@}
\section{cmpDisp}

@o cmpDisp.m
@{
function theRes=cmpDisp(rMat,sMat,tMat)
theRes=rMat-sMat*(tMat \ (sMat'));


@}
\section{cmpXHatEmpirical}

@o cmpXHatEmpirical.m
@{
function theRes=cmpXHatEmpirical(fs,gt,lilU,uu,yObs)
uInv=inv(uu);
preCalc=(gt')*uInv;
xZeroEst=(preCalc *gt)\ (preCalc *yObs);
theRes=(fs- lilU'*uInv*gt)*xZeroEst+ lilU'*uInv*yObs;

@}
\section{schurSoln}
@o schurSoln.m
@{

function xx=schurSoln(aa,bb)
theDim=size(aa,1);
theCols=size(bb,2);
[qqq,ttt]=schur(full(aa));
[qq,tt]=ordschur(qqq,ttt,ordeig(ttt)>1e-8); 
cc=qq'*bb;
zapThese=min(find(abs(diag(tt))<1e-10));
if(norm(cc(zapThese:end,:))>1e-10)
error('schurSoln: trouble zapping part of vector');
end
yy=zeros(theDim,theCols);
yy(1:(zapThese-1),:)=tt(1:(zapThese-1),1:(zapThese-1))\cc(1:(zapThese-1),:);
xx=qq*yy;
@}


\section{cmpBigC}

@o cmpBigC.m
@{
function theRes=cmpBigC(fs,gt,lilU,uInv)
theRes=schurSoln(((gt') * uInv * gt) , (fs' - (gt')*uInv*lilU));
@}
\section{cmpLStar}

@o cmpLStar.m
@{
function theRes=cmpLStar(fs,gt,lilU,uInv,bigC)
theRes=uInv*(lilU + gt * bigC);

@}
\section{cmpXHatNoInfo}

@o cmpXHatNoInfo.m
@{
function theRes=cmpXHatNoInfo(fs,gt,lilU,uInv,yObs,bigC)
%bigC=cmpBigC(fs,gt,lilU,uu);
lStar=cmpLStar(fs,gt,lilU,uInv,bigC);
theRes=(lStar')*yObs;

@}
\section{cmpDispNoInfo}

@o cmpDispNoInfo.m
@{
function [theRes1,theRes2]=cmpDispNoInfo(fs,gt,lilU,uInv,vss,bigC)
theRes1=[];%vss - lilU' * uInv * lilU;
theRes2=(vss - lilU' * uInv * lilU)+ bigC' * gt' * uInv * gt * bigC;

@}
\section{cmpLilU}

@o cmpLilU.m
@{
function theRes=cmpLilU(ff,hh,vv,tt,ss)
numRows=nCols(hh);
numCols=nRows(hh);
theRes=zeros(numCols*tt,numRows);
for ii=1:tt
theRes((1:numCols)+(ii-1)*numCols,:)=...
(cmpEDelCov(ff,hh,vv,ss,ii))';
end
@}
\section{cmpBigU}

@o cmpBigU.m
@{
function theRes=cmpBigU(ff,hh,vv,ww,tt)
hDim=nRows(hh);
theRes=zeros(tt*hDim);
vtt=[];%cmpECov(ff,vv,tt); never used i think
[newvss,newbigu,theWs]=cmpBigECov(ff,vv,hh,ww,tt);
for ii=1:tt
for jj=ii:tt
toGo=newCmpDelCov(ff,hh,vv,ww,vtt,ii,jj);
theRes((1:hDim)+(ii-1)*hDim,(1:hDim)+(jj-1)*hDim)=toGo;
if(jj>ii)
theRes((1:hDim)+(jj-1)*hDim,(1:hDim)+(ii-1)*hDim)=toGo';
end
end
end
@}
\section{prttnBigU}

@o prttnBigU.m
@{
function [theRes11,theRes12,theRes21,theRes22]=prttnBigU(uu,ss,tt)
blksize=nRows(uu)/tt;
oneRange=1:blksize*ss;
twoRange=((blksize*ss)+1):blksize*tt;
theRes11=uu(oneRange,oneRange);
theRes12=uu(oneRange,twoRange);
theRes21=uu(twoRange,oneRange);
theRes22=uu(twoRange,twoRange);
@}
\section{prttnLilU}

@o prttnLilU.m
@{
function [theRes11,theRes21]=prttnLilU(lilU,ss,tt)
blksize=nRows(lilU)/tt;
oneRange=1:(blksize*ss);
twoRange=((blksize*ss)+1):(blksize*tt);
theRes11=lilU(oneRange,:);
theRes21=lilU(twoRange,:);
@}
\section{prttnG}

@o prttnG.m
@{
function [theRes11,theRes21]=prttnG(gg,ss,tt)
blksize=nRows(gg)/tt;
oneRange=1:blksize*ss;
twoRange=((blksize*ss)+1):blksize*tt;
theRes11=gg(oneRange,:);
theRes21=gg(twoRange,:);

@}
\section{cmpRecursive}

@o cmpRecursive.m
@{
function [theRes1,theRes2,theRes3,theRes4,theRes5]=...
cmpRecursive(ff,hh,vv,ww,tt,ss)
fs=ff^ss;
gt=cmpGtt(ff,hh,tt);
lilU=cmpLilU(ff,hh,vv,tt,ss);
uu=cmpBigU(ff,hh,vv,ww,tt);
[uPrts0101,uPrts0102,uPrts0201,uPrts0202]=prttnBigU(uu,ss,tt);
[lilUPrts01,lilUPrts02]=prttnLilU(lilU,ss,tt);
[gPrts01,gPrts02]=prttnG(gt,ss,tt);
lStar1=cmpLStar(fs,gPrts01,lilUPrts01,uPrts0101);
lStar2=cmpLStar(gPrts02,gPrts01,uPrts0102,uPrts0101);
vss=cmpECov(ff,vv,ss);
theRes1=lStar1;
theRes2=lStar2;
theRes3=cmpCAAorDD(vss,lStar1,uPrts0101,lilUPrts01);
theRes4=cmpCAAorDD(uPrts0202,lStar2,uPrts0101,uPrts0102);
theRes5=cmpCAD(lStar1,lStar2,uPrts0101,uPrts0102,lilUPrts01,lilUPrts02);
@}
\section{cmpQRRecursive}

@o cmpQRRecursive.m
@{
function [theRes1,theRes2]=cmpQRRecursive(ff,hh,tt,ss)
fs=ff^ss;
gt=cmpGtt(ff,hh,tt);
[qq,rr,ee]=qr(gt);
theRes1=qq;
theRes2=ee*( rr \ qq);
@}
\section{cmpCAAorDD}

@o cmpCAAorDD.m
@{
function theRes=cmpCAAorDD(vss,lStar1,u11,lilU1)
theRes=(vss + lStar1' * u11 * lStar1 - lilU1' * lStar1 - lStar1' * lilU1);
@}
\section{cmpCAD}

@o cmpCAD.m
@{
function theRes=cmpCAD(lStar1,lStar2,u11,u12,lilU1,lilU2)
theRes=(lilU2' - lStar1' * u12 - lilU1' * lStar2 + lStar1' *  u11 * lStar2);
@}
\section{applyFilter}

@o applyFilter.m
@{
function [theRes1, theRes2]=applyFilter(lStar1,lStar2,caa,cdd,cad,yObs)
y1Range=(1:(nRows(lStar1)));
y2Range=nRows(lStar1)+(1:(nCols(lStar2)));
cddInv=inv(cdd);
theRes1=lStar1' * yObs(y1Range) +...
cad * cddInv * ...
(yObs(y2Range) - ...
lStar2' * yObs(y1Range));
theRes2=caa-cad * cddInv * cad';

@}
\section{toCompare}

@o toCompare.m
@{
function [cqr01,cqr02,theProd,recX,recDisp,...
xhat,disp,xhatnoinfo,xhatempirical,dispnoinfo]=...
toCompare(anX0,aPMat,ff,hh,vv,ww,bigS,bigT,varargin)
if(nargin>5)
yObs=varargin{1};
else
yObs=rand(bigT,1);
end
[cqr01,cqr02]=cmpQRRecursive(ff,hh,bigT,bigS);
theProd=cqr02 *yObs(1:(nRows(ww)*bigT));
[lStar1,lStar2,caa,cdd,cad]=cmpRecursive(ff,hh,vv,ww,bigT,bigS);
[recX,recDisp]=applyFilter(lStar1,lStar2,caa,cdd,cad,yObs);
[xhat,disp,xhatnoinfo,xhatempirical,dispnoinfo]=priorMeanVar(anX0,aPMat,ff,hh,vv,ww,bigT,bigS,yObs);
@}
\section{uncondCov}

@o uncondCov.m
@{

function varCov=uncondCov(theMat,shockCov)
%varCov=uncondCov(theMat,shockCov)
[srows,scols]=size(theMat);
varCov=reshape((speye((srows)^2)-kron(theMat,theMat)) \ ...
reshape(shockCov,(srows)^2,1),srows,srows);



@}
\section{simPath}

@o simPath.m
@{


function [xiVecs,etaVecs,xVecs,yVecs]=simPath(xZero,ff,hh,vv,ww,tt)
pp=nRows(ff);
qq=nRows(hh);
xiVecs=multNormDraw(vv,tt);
etaVecs=multNormDraw(ww,tt);
xVecs=zeros(pp,tt+1);
yVecs=zeros(qq,tt+1);
xVecs(:,1)=xZero;
yVecs(:,1)=hh * xZero+multNormDraw(ww,1);
for ii=1:tt
xVecs(:,ii+1)=ff *xVecs(:,ii)+xiVecs(:,ii);
yVecs(:,ii+1)=hh *xVecs(:,ii+1)+etaVecs(:,ii);
end
@}


\section{multNormDraw}

@o multNormDraw.m
@{
function theRes=multNormDraw(covMat,tt)
theDim=nRows(covMat);
theR=chol(covMat);
raw=randn(theDim,tt);
theRes=theR*raw;
@}
@o tryFilt.m
@{
function [xiVecs,etaVecs,xVecs,yVecs,mLogL,lhtsmooth,loglh,xLike]=...
tryFilt(theParams,anH,anX0,aV0,aW0,bigT,delS,delPS,trainvec)
%xpreds blue xhatKF green xhatSmth red xVecs black
%function [xiVecs,etaVecs,xVecs,yVecs,mLogL,lhtsmooth,loglh,xLike]=...
%tryFilt([a11 a12 a21 a22],[h0,h1],[x01 x02],[xvar1,cov;cov,xv2],[vary],bigT,...
%[x0guess0 ;x0guess1], [prec00 prec01;prec10 prec11])
nn=2;
anF=[reshape(theParams(1:4),nn,nn)];


[xiVecs,etaVecs,xVecs,yVecs]=simPath(anX0,anF,anH,aV0,aW0,bigT);



shat=zeros(nn,1)+delS;
xLike=gaussian_prob(xiVecs,zeros(size(shat,1),1),aV0,true);
xLike=sum(xLike);
pshat=uncondCov(anF,aV0)+delPS;

@}
@o tryFilt.m
@{
numRows=nRows(anX0);
xhatKF=zeros(numRows,bigT);
xhatSmth=zeros(numRows,bigT);
notShat=shat;
notPshat=pshat;
for ii=1:bigT;
%    [xhatKF(:,ii),pshat,lht(ii)]=notKfSmooth(yVecs(:,ii),anH,...
%notShat,pshat,anF,aV0);
%    [xhatKF(:,ii),notPshat,lht(ii)]=kf(yVecs(:,ii),anH,...
%notShat,notPshat,anF,aV0);
    [xhatKF(:,ii),notPshat,lht(ii)]=notKf(yVecs(:,ii),anH,...
notShat,notPshat,anF,aV0);
notShat=xhatKF(:,ii);
end;

lht=lht(trainvec(1):trainvec(2));
mLogL=-(bigT*nn*0.5)*log(2*pi)+sum(lht);


  [xhatSmth,pshatsmooth,lhtsmooth]=notKfSmooth(yVecs(:,1:bigT),...
anH,shat,pshat,anF,aV0);


@}
@o tryFilt.m
@{

[xpreds,ypreds,disp,theErrs,loglh]=cmpXHatsYHats(anF,anH,aV0,aW0,bigT,...
vec(yVecs));
loglh=sum(loglh);

theErrs=cmpPredErrs(anF,anH,aV0,aW0,bigT,vec(yVecs));
theRange=(1:size(xpreds,2));
figure(1)
plot(...
theRange,xpreds(2,theRange),'b',...
theRange,xhatKF(2,theRange),'g',...
theRange,xhatSmth(2,theRange),'r',...
theRange,xVecs(2,theRange),'k'...
)
figure(2)
plot(...
theRange,xpreds(1,theRange),'b',...
theRange,xhatKF(1,theRange),'g',...
theRange,xhatSmth(1,theRange),'r',...
theRange,xVecs(1,theRange),'k'...
)
figure(3)
plot(theErrs)
tilefig([1,2,3]);
@}
@o tryBigger.m
@{
function [mLogL,lhtsmooth,loglh,xLike]=tryBigger(theN,pathLen)
%[mLogL,lhtsmooth,loglh,xLike]=tryBigger(theN,pathLen)
mLogL=zeros(theN,1);
lhtsmooth=zeros(theN,1);
loglh=zeros(theN,1);
xLike=zeros(theN,1);
for ii=1:theN
[xiVecs,etaVecs,xVecs,yObs,mLogL(ii),lhtsmooth(ii),loglh(ii),xLike(ii)]=...
tryFilt([0.9,0.1 0.2,0.7],eye(2),[0;0],0.001*eye(2),0.001*eye(2),pathLen,...
zeros(2,1),zeros(2,2),[1 pathLen]);
end

corrcoef([mLogL,lhtsmooth,loglh,xLike])
figure(5)
plot(mLogL,xLike,'o');
figure(6)
plot(mLogL,xLike,'o');
figure(7)
plot(mLogL,xLike,'o');
tilefig([5 6 7])
@}
@o tryNow.m
@{
function [mLogL,lhtsmooth,loglh,xLike]=tryNow(theN,pathLen)

mLogL=zeros(theN,1);
lhtsmooth=zeros(theN,1);
loglh=zeros(theN,1);
xLike=zeros(theN,1);
for ii=1:theN
[xiVecs,etaVecs,xVecs,yObs,mLogL(ii),lhtsmooth(ii),loglh(ii),xLike(ii)]=...
tryFilt([0.9,0.1 0.2,0.7],[1 2],[0;0],eye(2),0.00,pathLen,...
zeros(2,1),zeros(2,2),[1 pathLen]);
end
@}


\section{tryIt}

@o tryIt.m
@{
anF=[0.9,0.1;0.2,0.7];
anH=[1,2];
anX0=[0.0;0.0];
aP0=eye(2);
aV0=eye(2);
aW0=[0.0000001];
bigT=15;
bigS=13;
bigT=10;
bigS=8;
bigT=30;
bigS=25;

[xiVecs,etaVecs,xVecs,yVecs]=simPath(anX0,anF,anH,aV0,aW0,bigT);

yObs=yVecs(:);

%[cqr01,cqr02,theProd,recX,recDisp,...
%xhat,disp,xhatnoinfo,xhatempirical,dispnoinfo]=...
%toCompare(anX0,aP0,anF,anH,aV0,aW0,bigS,bigT,yObs);

[xpreds,ypreds,disp,theErrs,loglh]=cmpXHatsYHats(anF,anH,aV0,aW0,bigT,yObs);

theErrs=cmpPredErrs(anF,anH,aV0,aW0,bigT,yObs);
theRange=(1:size(xpreds,2));
figure(1)
plot(theRange,xpreds(2,theRange),'b',theRange,xVecs(2,theRange),'r')
figure(2)
plot(theRange,xpreds(1,theRange),'b',theRange,xVecs(1,theRange),'r')
figure(3)
plot(theErrs)
@}

\end{document}
