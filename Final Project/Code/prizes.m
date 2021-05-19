clear
load prizes_bycomp

for i=1:size(competitionid),
    comp=competitionid(i)
    prize=[prize_1(i), prize_2(i), prize_3(i)];
    prize=prize/sum(prize);
    aux=exist(sprintf('%d',comp),'dir');
    if aux==7,
        %cd (sprintf('%d',comp))
        save(sprintf('%d/%s_%02d.mat', comp,'prize', comp),'prize')
    end
    if aux==0,
        mkdir(sprintf('%d',comp))
        save(sprintf('%d/%s_%02d.mat', comp,'prize', comp),'prize')
    end
end