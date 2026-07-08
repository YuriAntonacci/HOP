%% create (and possibly save) PID antichains and Mobius inversion matrices for N=2 sources
clear; close all; clc
salva='y';
filename='Mobius2.mat';

ac{1}=1; act{1}='[1]';
ac{2}=2; act{2}='[2]';
ac{3}=[1 2]; act{3}='[1,2]'; 
ac{4}{1}=1; ac{4}{2}=2; act{4}='[1],[2]';

cg(1,:)='U1';
cg(2,:)='U2';
cg(3,:)='Sy';
cg(4,:)='Rd';

M=[1 0 -1 0; 0 1 -1 0; 0 0 1 0; -1 -1 1 1];

if salva=='y'
    PID.type='2_sources_PID';
    PID.Ns=2;
    PID.ac=ac;
    PID.act=act;
    PID.M=M;
    PID.cg=cg;
    save(filename,'PID')
end


