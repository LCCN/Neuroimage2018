%ITPC computation according to Lachaux et al 1999, this script assumes the
%input is an EEG dataset (sensors o source space).

function [nsub,itpc] = itpc_taskswitching(data)

%% Filtering
%   La técnica elegida es coger filtros de 2 en 2 Hz entre 2 y 46, y definir
% las bandas a partir de ahí, de la forma tradicional 
% delta [2-4 Hz] para evitar efectos de borde acusados en el filtrado
% theta [4-8 Hz], alpha [8-12 Hz], beta [12-30 Hz], gamma [30-46 Hz]

% Primero se cargan los parámetros de los filtros, que han sido definidos
% en el script defiltro.m

Nb = 5;  % number of frequency bands
[~,Nf,Nt] = size(data); % Numero de puntos, fuentes y triales
Np2 = 850; % Samples

%Data Filtering
% data is a matrix for each subject and condition (Nsamples x Nsources x Ntrials), 

Ns=99;
phi=single(zeros(Nb,Np2,Nf,Nt)); % Matrix to save the phases
ws=5; Nw=floor(Np2/ws); % window length ; number of windows 
phi2=single(zeros(Nb,ws*Nt,Nf,Nw)); % Concatenated phases, Nw will be the time
phi3=phi2; % Surrogates
itpc=single(zeros(Nb,Nf,Nf,Nw)); % Output
temp=single(zeros(Nf,Nf)); nsub=itpc; % for the surrogates


% Phases computation  
parfor i=1:Nb
    dat_filt=data;
    % Remove half second before and after due to edge effect
    phi(i,:,:,:)=unwrap(angle(dat_filt(:,:,:))); 
end

%%  ITPC computation
   
clear dat_filt 

for i=1:Nw % For each window
    for j=1:Nt
         %Windows for each frequency bands, sources and trials
        phi2(:,(j-1)*ws+1:j*ws,:,i)=phi(:,(i-1)*ws+1:i*ws,:,j);
    end
end        

% ITPC
tic
for i=1:Nb
    a=phi2(i,:,:,:);  % Phases for each frequency band
    parfor j=1:Nw
        c=exp(1i*bsxfun(@minus,a(:,:,:,j),permute(a(:,:,:,j),[3,2,1,4])));
        itpc(i,:,:,j)=squeeze(abs(mean(c,2)));
    end
end
toc

%% Surrogates
orden=1:Nt; % sorted vector
count=0;

while count<Ns,
    shu=randperm(Nt); match=randperm(Nt)-orden; % Permutation and comparison with the original
    
    if sum(match==0)<Nt/4,  % Only a fourth of the trials can coincide in the permutation
        tic,count=count+1 %#ok<NOPRT>
     
        for i=1:Nw % For each window
            for j=1:Nt
                % Window for each  band, sources and random trials
                phi3(:,(shu(j)-1)*ws+1:shu(j)*ws,:,i)=phi(:,(i-1)*ws+1:i*ws,:,j);
            end
        end
        
        for i=1:Nb
            a=phi2(i,:,:,:); b=phi3(i,:,:,:); % Phases for each freq band
            
            parfor j=1:Nw
                c=exp(1i*bsxfun(@minus,a(:,:,:,j),permute(b(:,:,:,j),[3,2,1,4])));
                temp=squeeze(itpc(i,:,:,j))-squeeze(abs(mean(c,2))); % check which one is greater, the original or the surrogate      temp(temp>0)=0; temp(temp<0)=1; % 1 if the surrogate is greater 
                nsub(i,:,:,j)=squeeze(nsub(i,:,:,j))+temp; % Accumulate               
            end
            
        end      
        toc
    end
end
% Finally for 99 surrogates and 1 original, goes through with a p<0.05 the
% ones with an accumulate of 9 or less
return    
