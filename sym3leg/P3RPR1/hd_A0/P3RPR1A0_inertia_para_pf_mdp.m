% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [10x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RPR1G1P1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(10,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1P1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1P1A0_inertia_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1P1A0_inertia_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1P1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1P1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'P3RPR1G1P1A0_inertia_para_pf_mdp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:16
% EndTime: 2019-05-03 14:58:17
% DurationCPUTime: 1.11s
% Computational Cost: add. (1485->207), mult. (2308->381), div. (249->6), fcn. (1958->14), ass. (0->146)
t1041 = 2 * MDP(5);
t981 = legFrame(3,3);
t970 = sin(t981);
t973 = cos(t981);
t984 = sin(qJ(1,3));
t987 = cos(qJ(1,3));
t939 = t970 * t987 + t973 * t984;
t940 = -t970 * t984 + t973 * t987;
t1001 = koppelP(3,1);
t991 = xP(3);
t979 = sin(t991);
t980 = cos(t991);
t998 = koppelP(3,2);
t945 = -t1001 * t979 - t980 * t998;
t948 = t1001 * t980 - t979 * t998;
t992 = 0.1e1 / qJ(2,3);
t894 = (t939 * t948 + t940 * t945) * t992;
t1040 = pkin(1) * t894;
t982 = legFrame(2,3);
t971 = sin(t982);
t974 = cos(t982);
t985 = sin(qJ(1,2));
t988 = cos(qJ(1,2));
t941 = t971 * t988 + t974 * t985;
t942 = -t971 * t985 + t974 * t988;
t1002 = koppelP(2,1);
t999 = koppelP(2,2);
t946 = -t1002 * t979 - t980 * t999;
t949 = t1002 * t980 - t979 * t999;
t994 = 0.1e1 / qJ(2,2);
t895 = (t941 * t949 + t942 * t946) * t994;
t1039 = pkin(1) * t895;
t983 = legFrame(1,3);
t972 = sin(t983);
t975 = cos(t983);
t986 = sin(qJ(1,1));
t989 = cos(qJ(1,1));
t943 = t972 * t989 + t975 * t986;
t944 = -t972 * t986 + t975 * t989;
t1000 = koppelP(1,2);
t1003 = koppelP(1,1);
t947 = -t1000 * t980 - t1003 * t979;
t950 = -t1000 * t979 + t1003 * t980;
t996 = 0.1e1 / qJ(2,1);
t896 = (t943 * t950 + t944 * t947) * t996;
t1038 = pkin(1) * t896;
t1037 = pkin(1) * t939;
t1036 = pkin(1) * t940;
t1035 = pkin(1) * t941;
t1034 = pkin(1) * t942;
t1033 = pkin(1) * t943;
t1032 = pkin(1) * t944;
t952 = t1001 * t984 - t987 * t998;
t953 = t1001 * t987 + t984 * t998;
t891 = (t952 * t980 - t953 * t979) * t973 + t970 * (t952 * t979 + t953 * t980);
t1031 = t891 * t894;
t954 = t1002 * t985 - t988 * t999;
t955 = t1002 * t988 + t985 * t999;
t892 = (t954 * t980 - t955 * t979) * t974 + t971 * (t954 * t979 + t955 * t980);
t1030 = t892 * t895;
t956 = -t1000 * t989 + t1003 * t986;
t957 = t1000 * t986 + t1003 * t989;
t893 = (t956 * t980 - t957 * t979) * t975 + t972 * (t956 * t979 + t957 * t980);
t1029 = t893 * t896;
t990 = pkin(1) + pkin(2);
t958 = -qJ(2,3) * t987 + t984 * t990;
t961 = qJ(2,3) * t984 + t987 * t990;
t921 = t958 * t973 + t961 * t970;
t1028 = (-t921 + 0.2e1 * t1037) * t992 ^ 2;
t924 = -t958 * t970 + t961 * t973;
t905 = (-t924 + 0.2e1 * t1036) * t992;
t1027 = t905 * t992;
t959 = -qJ(2,2) * t988 + t985 * t990;
t962 = qJ(2,2) * t985 + t988 * t990;
t922 = t959 * t974 + t962 * t971;
t1026 = (-t922 + 0.2e1 * t1035) * t994 ^ 2;
t925 = -t959 * t971 + t962 * t974;
t909 = (-t925 + 0.2e1 * t1034) * t994;
t1025 = t909 * t994;
t960 = -qJ(2,1) * t989 + t986 * t990;
t963 = qJ(2,1) * t986 + t989 * t990;
t923 = t960 * t975 + t963 * t972;
t1024 = (-t923 + 0.2e1 * t1033) * t996 ^ 2;
t926 = -t960 * t972 + t963 * t975;
t913 = (-t926 + 0.2e1 * t1032) * t996;
t1023 = t913 * t996;
t1022 = t939 * t992;
t1007 = qJ(2,3) ^ 2;
t993 = 0.1e1 / t1007;
t1021 = t939 * t993;
t1020 = t940 * t992;
t1019 = t940 * t993;
t1018 = t941 * t994;
t1006 = qJ(2,2) ^ 2;
t995 = 0.1e1 / t1006;
t1017 = t941 * t995;
t1016 = t942 * t994;
t1015 = t942 * t995;
t1014 = t943 * t996;
t1005 = qJ(2,1) ^ 2;
t997 = 0.1e1 / t1005;
t1013 = t943 * t997;
t1012 = t944 * t996;
t1011 = t944 * t997;
t1010 = (t921 * t948 + t924 * t945) * t992;
t1009 = (t922 * t949 + t925 * t946) * t994;
t1008 = (t923 * t950 + t926 * t947) * t996;
t1004 = pkin(1) ^ 2;
t978 = t1004 + t1005;
t977 = t1004 + t1006;
t976 = t1004 + t1007;
t969 = qJ(2,1) * t1003 + t1000 * t990;
t968 = qJ(2,2) * t1002 + t990 * t999;
t967 = qJ(2,3) * t1001 + t990 * t998;
t966 = -qJ(2,1) * t1000 + t1003 * t990;
t965 = -qJ(2,2) * t999 + t1002 * t990;
t964 = -qJ(2,3) * t998 + t1001 * t990;
t951 = (t979 ^ 2 + t980 ^ 2) * MDP(10);
t938 = t944 ^ 2;
t937 = t943 ^ 2;
t936 = t942 ^ 2;
t935 = t941 ^ 2;
t934 = t940 ^ 2;
t933 = t939 ^ 2;
t932 = t966 * t986 - t969 * t989;
t931 = t965 * t985 - t968 * t988;
t930 = t964 * t984 - t967 * t987;
t929 = t966 * t989 + t969 * t986;
t928 = t965 * t988 + t968 * t985;
t927 = t964 * t987 + t967 * t984;
t914 = (t926 - t1032) * t996;
t912 = (t923 - t1033) * t996;
t910 = (t925 - t1034) * t994;
t908 = (t922 - t1035) * t994;
t906 = (t924 - t1036) * t992;
t904 = (t921 - t1037) * t992;
t902 = (-pkin(1) * t926 + t944 * t978) * t996;
t901 = (-pkin(1) * t923 + t943 * t978) * t996;
t900 = (-pkin(1) * t925 + t942 * t977) * t994;
t899 = (-pkin(1) * t922 + t941 * t977) * t994;
t898 = (-pkin(1) * t924 + t940 * t976) * t992;
t897 = (-pkin(1) * t921 + t939 * t976) * t992;
t890 = (-t929 * t979 + t932 * t980) * t975 + (t929 * t980 + t932 * t979) * t972;
t889 = (-t928 * t979 + t931 * t980) * t974 + (t928 * t980 + t931 * t979) * t971;
t888 = (-t927 * t979 + t930 * t980) * t973 + (t927 * t980 + t930 * t979) * t970;
t1 = [(t934 * t993 + t936 * t995 + t938 * t997) * MDP(1) + ((-t926 * t997 + t1023) * t944 + (-t925 * t995 + t1025) * t942 + (-t924 * t993 + t1027) * t940) * MDP(4) + (t934 * t992 + t936 * t994 + t938 * t996) * t1041 + ((t902 * t944 + t914 * t926) * t996 + (t900 * t942 + t910 * t925) * t994 + (t898 * t940 + t906 * t924) * t992) * MDP(6) + t951; (t1011 * t943 + t1015 * t941 + t1019 * t939) * MDP(1) + (-t1011 * t923 + t1014 * t913 - t1015 * t922 + t1018 * t909 - t1019 * t921 + t1022 * t905) * MDP(4) + (t1012 * t943 + t1016 * t941 + t1020 * t939) * t1041 + ((t902 * t943 + t914 * t923) * t996 + (t900 * t941 + t910 * t922) * t994 + (t898 * t939 + t906 * t921) * t992) * MDP(6); (t933 * t993 + t935 * t995 + t937 * t997) * MDP(1) + ((-t923 * t997 + t1024) * t943 + (-t922 * t995 + t1026) * t941 + (-t921 * t993 + t1028) * t939) * MDP(4) + (t933 * t992 + t935 * t994 + t937 * t996) * t1041 + ((t901 * t943 + t912 * t923) * t996 + (t899 * t941 + t908 * t922) * t994 + (t897 * t939 + t904 * t921) * t992) * MDP(6) + t951; (t1011 * t893 + t1015 * t892 + t1019 * t891) * MDP(1) + (-t1011 * t890 - t1015 * t889 - t1019 * t888 + t1023 * t893 + t1025 * t892 + t1027 * t891) * MDP(4) + (t1012 * t893 + t1016 * t892 + t1020 * t891) * t1041 + ((t890 * t914 + t893 * t902) * t996 + (t889 * t910 + t892 * t900) * t994 + (t888 * t906 + t891 * t898) * t992) * MDP(6) - t979 * MDP(8) - t980 * MDP(9); (t1013 * t893 + t1017 * t892 + t1021 * t891) * MDP(1) + (-t1013 * t890 - t1017 * t889 - t1021 * t888 + t1024 * t893 + t1026 * t892 + t1028 * t891) * MDP(4) + (t1014 * t893 + t1018 * t892 + t1022 * t891) * t1041 + ((t890 * t912 + t893 * t901) * t996 + (t889 * t908 + t892 * t899) * t994 + (t888 * t904 + t891 * t897) * t992) * MDP(6) + t980 * MDP(8) - t979 * MDP(9); (t1029 + t1030 + t1031) * t1041 + MDP(7) + (MDP(1) * t1029 + (t893 * (-t1008 + 0.2e1 * t1038) - t890 * t896) * MDP(4) + (t893 * (-pkin(1) * t1008 + t978 * t896) + t890 * (t1008 - t1038)) * MDP(6)) * t996 + (MDP(1) * t1030 + (t892 * (-t1009 + 0.2e1 * t1039) - t889 * t895) * MDP(4) + (t892 * (-pkin(1) * t1009 + t977 * t895) + t889 * (t1009 - t1039)) * MDP(6)) * t994 + (MDP(1) * t1031 + (t891 * (-t1010 + 0.2e1 * t1040) - t888 * t894) * MDP(4) + (t891 * (-pkin(1) * t1010 + t976 * t894) + t888 * (t1010 - t1040)) * MDP(6)) * t992;];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1), t1(2), t1(4); t1(2), t1(3), t1(5); t1(4), t1(5), t1(6);];
MMX  = res;
