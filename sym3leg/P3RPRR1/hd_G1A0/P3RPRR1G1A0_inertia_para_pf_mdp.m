% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRR1G1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RPRR1G1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1A0_inertia_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3RPRR1G1A0_inertia_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:21
% EndTime: 2020-03-09 21:23:22
% DurationCPUTime: 0.79s
% Computational Cost: add. (1727->140), mult. (1113->238), div. (276->7), fcn. (984->29), ass. (0->126)
t1025 = MDP(4) * pkin(1) ^ 2 + MDP(1);
t963 = sin(pkin(7));
t1024 = pkin(1) * t963;
t960 = legFrame(3,3) + qJ(1,3);
t956 = pkin(7) + t960;
t953 = qJ(3,3) + t956;
t947 = sin(t953);
t911 = -pkin(1) * sin(t960) - pkin(2) * sin(t956) - pkin(3) * t947;
t968 = sin(qJ(3,3));
t938 = pkin(1) * sin(pkin(7) + qJ(3,3)) + t968 * pkin(2);
t932 = 0.1e1 / t938;
t1022 = t911 * t932;
t961 = legFrame(2,3) + qJ(1,2);
t957 = pkin(7) + t961;
t954 = qJ(3,2) + t957;
t948 = sin(t954);
t912 = -pkin(1) * sin(t961) - pkin(2) * sin(t957) - pkin(3) * t948;
t969 = sin(qJ(3,2));
t939 = pkin(1) * sin(pkin(7) + qJ(3,2)) + t969 * pkin(2);
t934 = 0.1e1 / t939;
t1021 = t912 * t934;
t962 = legFrame(1,3) + qJ(1,1);
t958 = pkin(7) + t962;
t955 = qJ(3,1) + t958;
t949 = sin(t955);
t913 = -pkin(1) * sin(t962) - pkin(2) * sin(t958) - pkin(3) * t949;
t970 = sin(qJ(3,1));
t940 = pkin(1) * sin(pkin(7) + qJ(3,1)) + t970 * pkin(2);
t936 = 0.1e1 / t940;
t1020 = t913 * t936;
t950 = cos(t953);
t914 = -pkin(1) * cos(t960) - pkin(2) * cos(t956) - pkin(3) * t950;
t1019 = t914 * t932;
t951 = cos(t954);
t915 = -pkin(1) * cos(t961) - pkin(2) * cos(t957) - pkin(3) * t951;
t1018 = t915 * t934;
t952 = cos(t955);
t916 = -pkin(1) * cos(t962) - pkin(2) * cos(t958) - pkin(3) * t952;
t1017 = t916 * t936;
t974 = 0.1e1 / pkin(3);
t1016 = t932 * t974;
t933 = 0.1e1 / t938 ^ 2;
t1015 = t933 * t947;
t1014 = t933 * t950;
t1013 = t934 * t974;
t935 = 0.1e1 / t939 ^ 2;
t1012 = t935 * t948;
t1011 = t935 * t951;
t1010 = t936 * t974;
t937 = 0.1e1 / t940 ^ 2;
t1009 = t937 * t949;
t1008 = t937 * t952;
t917 = t947 * t932;
t918 = t948 * t934;
t919 = t949 * t936;
t920 = t950 * t932;
t921 = t951 * t934;
t922 = t952 * t936;
t1007 = t963 * t968;
t1006 = t963 * t969;
t1005 = t963 * t970;
t1004 = t1025 * (t949 * t1008 + t948 * t1011 + t947 * t1014);
t1003 = 2 * MDP(7);
t964 = cos(pkin(7));
t971 = cos(qJ(3,3));
t923 = -t971 * pkin(2) + (-t964 * t971 + t1007) * pkin(1);
t1002 = t923 * t917;
t1001 = t923 * t920;
t972 = cos(qJ(3,2));
t924 = -t972 * pkin(2) + (-t964 * t972 + t1006) * pkin(1);
t1000 = t924 * t918;
t999 = t924 * t921;
t973 = cos(qJ(3,1));
t925 = -t973 * pkin(2) + (-t964 * t973 + t1005) * pkin(1);
t998 = t925 * t919;
t997 = t925 * t922;
t959 = t964 * pkin(1) + pkin(2);
t926 = -pkin(1) * t1007 + t971 * t959;
t996 = t926 * t1015;
t995 = t926 * t1014;
t927 = -pkin(1) * t1006 + t972 * t959;
t994 = t927 * t1012;
t993 = t927 * t1011;
t928 = -pkin(1) * t1005 + t973 * t959;
t992 = t928 * t1009;
t991 = t928 * t1008;
t929 = t971 * t1024 + t968 * t959;
t990 = t929 * t917;
t989 = t929 * t920;
t988 = t929 * t1015;
t987 = t929 * t1014;
t930 = t972 * t1024 + t969 * t959;
t986 = t930 * t918;
t985 = t930 * t921;
t984 = t930 * t1012;
t983 = t930 * t1011;
t931 = t973 * t1024 + t970 * t959;
t982 = t931 * t919;
t981 = t931 * t922;
t980 = t931 * t1009;
t979 = t931 * t1008;
t910 = t916 * t1010;
t909 = t915 * t1013;
t908 = t914 * t1016;
t907 = t913 * t1010;
t906 = t912 * t1013;
t905 = t911 * t1016;
t904 = t922 + t910;
t903 = t921 + t909;
t902 = t920 + t908;
t901 = t919 + t907;
t900 = t918 + t906;
t899 = t917 + t905;
t898 = t910 + 0.2e1 * t922;
t897 = t909 + 0.2e1 * t921;
t896 = t908 + 0.2e1 * t920;
t895 = t907 + 0.2e1 * t919;
t894 = t906 + 0.2e1 * t918;
t893 = t905 + 0.2e1 * t917;
t892 = t922 + t910 / 0.2e1;
t891 = t921 + t909 / 0.2e1;
t890 = t920 + t908 / 0.2e1;
t889 = t919 + t907 / 0.2e1;
t888 = t918 + t906 / 0.2e1;
t887 = t917 + t905 / 0.2e1;
t1 = [(t902 * t920 + t903 * t921 + t904 * t922) * MDP(5) + (-t896 * t1001 - t897 * t999 - t898 * t997) * MDP(6) + (-t890 * t989 - t891 * t985 - t892 * t981) * t1003 + MDP(8) + ((t904 * t1017 + t903 * t1018 + t902 * t1019) * MDP(5) + (t914 * t995 + t915 * t993 + t916 * t991) * MDP(6) + (-t914 * t987 - t915 * t983 - t916 * t979) * MDP(7)) * t974 + t1025 * (t933 * t950 ^ 2 + t935 * t951 ^ 2 + t937 * t952 ^ 2); (t899 * t920 + t900 * t921 + t901 * t922) * MDP(5) + (-t893 * t1001 - t894 * t999 - t895 * t997) * MDP(6) + (-t887 * t989 - t888 * t985 - t889 * t981) * t1003 + ((t901 * t1017 + t900 * t1018 + t899 * t1019) * MDP(5) + (t914 * t996 + t915 * t994 + t916 * t992) * MDP(6) + (-t914 * t988 - t915 * t984 - t916 * t980) * MDP(7)) * t974 + t1004; 0; (t902 * t917 + t903 * t918 + t904 * t919) * MDP(5) + (-t897 * t1000 - t896 * t1002 - t898 * t998) * MDP(6) + (-t890 * t990 - t891 * t986 - t892 * t982) * t1003 + ((t904 * t1020 + t903 * t1021 + t902 * t1022) * MDP(5) + (t911 * t995 + t912 * t993 + t913 * t991) * MDP(6) + (-t911 * t987 - t912 * t983 - t913 * t979) * MDP(7)) * t974 + t1004; (t899 * t917 + t900 * t918 + t901 * t919) * MDP(5) + (-t894 * t1000 - t893 * t1002 - t895 * t998) * MDP(6) + (-t887 * t990 - t888 * t986 - t889 * t982) * t1003 + MDP(8) + ((t901 * t1020 + t900 * t1021 + t899 * t1022) * MDP(5) + (t911 * t996 + t912 * t994 + t913 * t992) * MDP(6) + (-t911 * t988 - t912 * t984 - t913 * t980) * MDP(7)) * t974 + t1025 * (t933 * t947 ^ 2 + t935 * t948 ^ 2 + t937 * t949 ^ 2); 0; 0; 0; (3 * MDP(4)) + MDP(8);];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
