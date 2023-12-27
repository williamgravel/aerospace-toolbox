function [CA,G_1,G_2i,G_2] = generateCACode(PRN,len,opts)
arguments
    PRN (1,1) double
    len (1,1) double
    opts.CodeType (1,1) string {mustBeMember(opts.CodeType,["binary","signed"])} = "binary"
end

sv = selectSVCode(PRN);

G_1 = generateGCode(len,10,[3,10]);
[G_2, G_2i] = generateGCode(len,10,[2,3,6,8,9,10],sv{1,{'PS 1','PS 2'}});

CA = bitxor(G_1,G_2i);

if strcmp(opts.CodeType,"signed")
    CA = -2*CA + 1;
    G_1 = -2*G_1 + 1;
    G_2i = -2*G_2i + 1;
    G_2 = -2*G_2 + 1;
end

end
