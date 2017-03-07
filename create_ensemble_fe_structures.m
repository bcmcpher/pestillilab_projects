%% create ensemble fe structures
% Brent McPherson
%

% create streamline indices
fibs = 1:500000;

% randomly sample streamlines
rfib = randsample(fibs, 500000);

% either pull 100k from each parameter - within PROB/DETR
x1 = rfib(1:100000);
x2 = rfib(100001:200000);
x3 = rfib(200001:300000);
x4 = rfib(300001:400000);
x5 = rfib(400001:500000);
% five total ensemble lmax with 600k, 1.2, or 1.3 million streamlines

% pull 50k from each parameter - across PROB/DETR
x01 = rfib(1:50000);
x02 = rfib(50001:100000);
x03 = rfib(100001:150000);
x04 = rfib(150001:200000);
x05 = rfib(200001:250000);
x06 = rfib(250001:300000);
x07 = rfib(300001:350000);
x08 = rfib(350001:400000);
x09 = rfib(400001:450000);
x10 = rfib(450001:500000);
% 10 total ensemble lmax with 300k, 600k, or 650k

% or randomly sample 500k equally from 13 sources
% 7T data only has 9...

% 38462 is the closest fraction of 13 to get to equal parts
y01 = rfib(1:38463);
y02 = rfib(38464:76925);
y03 = rfib(76926:115387);
y04 = rfib(115388:153849);
y05 = rfib(153850:192311);
y06 = rfib(192312:230773);
y07 = rfib(230774:269235);
y08 = rfib(269236:307697);
y09 = rfib(307698:346159);
y10 = rfib(346160:384261);
y11 = rfib(384262:423083);
y12 = rfib(423084:461545);
y13 = rfib(461546:500000);

% 55556 is the closest fraction for 9 equal parts
z01 = rfib(1:55556);
z02 = rfib(55557:111113);
z03 = rfib(111114:166669);
z04 = rfib(166670:222225);
z05 = rfib(222226:277781);
z06 = rfib(277782:333337);
z07 = rfib(333338:388893);
z08 = rfib(388894:444449);
z09 = rfib(444450:500000);

