% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRR1G1P1A0
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3PRRR1G1P1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1P1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1P1A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1P1A0_inertia_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1P1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1P1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR1G1P1A0_inertia_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:57
% EndTime: 2020-03-09 21:14:58
% DurationCPUTime: 0.79s
% Computational Cost: add. (3059->120), mult. (1469->236), div. (456->6), fcn. (2112->24), ass. (0->123)
t934 = pkin(7) + qJ(2,1);
t962 = qJ(3,1) + t934;
t922 = cos(t962);
t925 = sin(t934);
t928 = cos(t934);
t947 = sin(t962);
t1016 = 0.1e1 / (-t925 * t922 + t947 * t928);
t1015 = legFrame(1,3);
t931 = cos(t1015);
t965 = sin(t1015);
t915 = t965 * t922 + t931 * t947;
t991 = t1016 * t915;
t933 = pkin(7) + qJ(2,2);
t961 = qJ(3,2) + t933;
t921 = cos(t961);
t924 = sin(t933);
t927 = cos(t933);
t946 = sin(t961);
t1017 = 0.1e1 / (-t924 * t921 + t946 * t927);
t1014 = legFrame(2,3);
t930 = cos(t1014);
t964 = sin(t1014);
t913 = t964 * t921 + t930 * t946;
t995 = t1017 * t913;
t932 = pkin(7) + qJ(2,3);
t960 = qJ(3,3) + t932;
t920 = cos(t960);
t923 = sin(t932);
t926 = cos(t932);
t945 = sin(t960);
t1018 = 0.1e1 / (-t923 * t920 + t945 * t926);
t1013 = legFrame(3,3);
t929 = cos(t1013);
t963 = sin(t1013);
t911 = t963 * t920 + t929 * t945;
t999 = t1018 * t911;
t1012 = MDP(2) / pkin(2) ^ 2;
t896 = pkin(2) * (t923 * t929 + t963 * t926) + t911 * pkin(3);
t1011 = t896 * t1018;
t897 = pkin(2) * (t924 * t930 + t964 * t927) + t913 * pkin(3);
t1010 = t897 * t1017;
t898 = pkin(2) * (t925 * t931 + t965 * t928) + t915 * pkin(3);
t1009 = t898 * t1016;
t912 = t929 * t920 - t945 * t963;
t899 = -pkin(2) * (t923 * t963 - t929 * t926) + t912 * pkin(3);
t1008 = t899 * t1018;
t914 = t930 * t921 - t946 * t964;
t900 = -pkin(2) * (t924 * t964 - t930 * t927) + t914 * pkin(3);
t1007 = t900 * t1017;
t916 = t931 * t922 - t947 * t965;
t901 = -pkin(2) * (t925 * t965 - t931 * t928) + t916 * pkin(3);
t1006 = t901 * t1016;
t1005 = t1018 ^ 2;
t942 = 0.1e1 / pkin(2);
t1004 = t1018 * t942;
t1003 = t1017 ^ 2;
t1002 = t1017 * t942;
t1001 = t1016 ^ 2;
t1000 = t1016 * t942;
t998 = t1018 * t912;
t997 = t1018 * sin(qJ(3,3));
t996 = t1018 * cos(qJ(3,3));
t994 = t1017 * t914;
t993 = t1017 * sin(qJ(3,2));
t992 = t1017 * cos(qJ(3,2));
t990 = t1016 * t916;
t989 = t1016 * sin(qJ(3,1));
t988 = t1016 * cos(qJ(3,1));
t941 = 0.1e1 / pkin(3);
t987 = t941 * t942;
t986 = t911 * t997;
t985 = t911 * t996;
t984 = t912 * t997;
t983 = t912 * t996;
t982 = t1018 * t987;
t981 = t913 * t993;
t980 = t913 * t992;
t979 = t914 * t993;
t978 = t914 * t992;
t977 = t1017 * t987;
t976 = t915 * t989;
t975 = t915 * t988;
t974 = t916 * t989;
t973 = t916 * t988;
t972 = t1016 * t987;
t971 = t911 * t1004;
t970 = t912 * t1004;
t969 = t913 * t1002;
t968 = t914 * t1002;
t967 = t915 * t1000;
t966 = t916 * t1000;
t959 = t1018 * t986;
t958 = t1018 * t985;
t957 = t1018 * t984;
t956 = t1018 * t983;
t955 = t1017 * t981;
t954 = t1017 * t980;
t953 = t1017 * t979;
t952 = t1017 * t978;
t951 = t1016 * t976;
t950 = t1016 * t975;
t949 = t1016 * t974;
t948 = t1016 * t973;
t944 = (t990 * t991 + t994 * t995 + t998 * t999) * t1012;
t895 = t901 * t972;
t894 = t900 * t977;
t893 = t899 * t982;
t892 = t898 * t972;
t891 = t897 * t977;
t890 = t896 * t982;
t889 = -t895 + t966;
t888 = -t892 + t967;
t887 = -t894 + t968;
t886 = -t891 + t969;
t885 = -t893 + t970;
t884 = -t890 + t971;
t883 = -t895 + 0.2e1 * t966;
t882 = -t894 + 0.2e1 * t968;
t881 = -t893 + 0.2e1 * t970;
t880 = -t892 + 0.2e1 * t967;
t879 = -t891 + 0.2e1 * t969;
t878 = -t890 + 0.2e1 * t971;
t1 = [(t881 * t983 + t882 * t978 + t883 * t973) * MDP(6) + (-t881 * t984 - t882 * t979 - t883 * t974) * MDP(7) + MDP(8) + (t916 ^ 2 * t1001 + t914 ^ 2 * t1003 + t912 ^ 2 * t1005) * t1012 + ((t885 * t998 + t887 * t994 + t889 * t990) * MDP(5) + ((-t889 * t1006 - t887 * t1007 - t885 * t1008) * MDP(5) + (-t899 * t956 - t900 * t952 - t901 * t948) * MDP(6) + (t899 * t957 + t900 * t953 + t901 * t949) * MDP(7)) * t941) * t942; (t878 * t983 + t879 * t978 + t880 * t973) * MDP(6) + (-t878 * t984 - t879 * t979 - t880 * t974) * MDP(7) + t944 + ((t884 * t998 + t886 * t994 + t888 * t990) * MDP(5) + ((-t888 * t1006 - t886 * t1007 - t884 * t1008) * MDP(5) + (-t899 * t958 - t900 * t954 - t901 * t950) * MDP(6) + (t899 * t959 + t900 * t955 + t901 * t951) * MDP(7)) * t941) * t942; 0; (t881 * t985 + t882 * t980 + t883 * t975) * MDP(6) + (-t881 * t986 - t882 * t981 - t883 * t976) * MDP(7) + t944 + ((t885 * t999 + t887 * t995 + t889 * t991) * MDP(5) + ((-t889 * t1009 - t887 * t1010 - t885 * t1011) * MDP(5) + (-t896 * t956 - t897 * t952 - t898 * t948) * MDP(6) + (t896 * t957 + t897 * t953 + t898 * t949) * MDP(7)) * t941) * t942; (t878 * t985 + t879 * t980 + t880 * t975) * MDP(6) + (-t878 * t986 - t879 * t981 - t880 * t976) * MDP(7) + MDP(8) + (t915 ^ 2 * t1001 + t913 ^ 2 * t1003 + t911 ^ 2 * t1005) * t1012 + ((t884 * t999 + t886 * t995 + t888 * t991) * MDP(5) + ((-t888 * t1009 - t886 * t1010 - t884 * t1011) * MDP(5) + (-t896 * t958 - t897 * t954 - t898 * t950) * MDP(6) + (t896 * t959 + t897 * t955 + t898 * t951) * MDP(7)) * t941) * t942; 0; 0; 0; (3 * MDP(1)) + MDP(8);];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
