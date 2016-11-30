

function sensregnames = estimate_sens_names(sensregstr)

    sensregnames = sensregstr;
    sensregnames = strrep(sensregnames,'.*','');
    sensregnames = strrep(sensregnames,'sqrt(','');
    sensregnames = strrep(sensregnames,')','');

    % Simplify naming conventions
    sensregnames = strrep(sensregnames,'sl1sr1','SAs');
    sensregnames = strrep(sensregnames,'sla1sra1','SAd');
    sensregnames = strrep(sensregnames,'sl2sr2','SBs');
    sensregnames = strrep(sensregnames,'sla2sra2','SBd');
    sensregnames = strrep(sensregnames,'sl9sr9','FRs');
    sensregnames = strrep(sensregnames,'sla9sra9','FRd');
end