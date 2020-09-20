% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRR1G2P2A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x8]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRR1G2P2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2P2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:25:07
% EndTime: 2020-03-09 21:25:08
% DurationCPUTime: 0.53s
% Computational Cost: add. (933->139), mult. (651->170), div. (81->4), fcn. (570->63), ass. (0->116)
t1030 = 0.1e1 / pkin(3);
t1072 = t1030 / 0.2e1;
t1022 = legFrame(1,2);
t1039 = qJ(1,1) + pkin(7);
t994 = t1022 + t1039;
t988 = qJ(3,1) + t994;
t995 = -t1022 + t1039;
t989 = qJ(3,1) + t995;
t962 = -cos(t989) + cos(t988);
t1036 = pkin(7) + qJ(3,1);
t971 = 0.1e1 / (pkin(1) * sin(t1036) + sin(qJ(3,1)) * pkin(2));
t1043 = t962 * t971;
t1071 = t1043 / 0.2e1;
t1021 = legFrame(2,2);
t1038 = qJ(1,2) + pkin(7);
t992 = t1021 + t1038;
t986 = qJ(3,2) + t992;
t993 = -t1021 + t1038;
t987 = qJ(3,2) + t993;
t961 = -cos(t987) + cos(t986);
t1035 = pkin(7) + qJ(3,2);
t970 = 0.1e1 / (pkin(1) * sin(t1035) + sin(qJ(3,2)) * pkin(2));
t1044 = t961 * t970;
t1070 = t1044 / 0.2e1;
t1020 = legFrame(3,2);
t1037 = qJ(1,3) + pkin(7);
t990 = t1020 + t1037;
t984 = qJ(3,3) + t990;
t991 = -t1020 + t1037;
t985 = qJ(3,3) + t991;
t960 = -cos(t985) + cos(t984);
t1034 = pkin(7) + qJ(3,3);
t969 = 0.1e1 / (pkin(1) * sin(t1034) + sin(qJ(3,3)) * pkin(2));
t1045 = t960 * t969;
t1069 = t1045 / 0.2e1;
t959 = sin(t988) + sin(t989);
t1046 = t959 * t971;
t1068 = t1046 / 0.2e1;
t958 = sin(t986) + sin(t987);
t1047 = t958 * t970;
t1067 = t1047 / 0.2e1;
t957 = sin(t984) + sin(t985);
t1048 = t957 * t969;
t1066 = t1048 / 0.2e1;
t1025 = sin(qJ(1,1));
t1028 = cos(qJ(1,1));
t1013 = sin(t1022);
t1016 = cos(t1022);
t968 = t1016 * g(1) - t1013 * g(2);
t1031 = g(3) * t1025 - t968 * t1028;
t1065 = t1031 * t1043;
t1064 = t1031 * t1046;
t1024 = sin(qJ(1,2));
t1027 = cos(qJ(1,2));
t1012 = sin(t1021);
t1015 = cos(t1021);
t967 = t1015 * g(1) - t1012 * g(2);
t1032 = g(3) * t1024 - t967 * t1027;
t1063 = t1032 * t1044;
t1062 = t1032 * t1047;
t1023 = sin(qJ(1,3));
t1026 = cos(qJ(1,3));
t1011 = sin(t1020);
t1014 = cos(t1020);
t966 = t1014 * g(1) - t1011 * g(2);
t1033 = g(3) * t1023 - t966 * t1026;
t1061 = t1033 * t1045;
t1060 = t1033 * t1048;
t1004 = qJ(1,1) + t1036;
t1001 = cos(t1004);
t1040 = t1001 * t971;
t1003 = qJ(1,2) + t1035;
t1000 = cos(t1003);
t1041 = t1000 * t970;
t1002 = qJ(1,3) + t1034;
t999 = cos(t1002);
t1042 = t969 * t999;
t1059 = t1031 * t1040 + t1032 * t1041 + t1033 * t1042;
t1058 = pkin(1) / 0.2e1;
t996 = sin(t1002);
t939 = g(3) * t996 - t966 * t999;
t1057 = t939 * t969;
t997 = sin(t1003);
t940 = g(3) * t997 - t967 * t1000;
t1056 = t940 * t970;
t998 = sin(t1004);
t941 = g(3) * t998 - t968 * t1001;
t1055 = t941 * t971;
t942 = g(3) * t999 + t966 * t996;
t1054 = t942 * t969;
t943 = g(3) * t1000 + t967 * t997;
t1053 = t943 * t970;
t944 = g(3) * t1001 + t968 * t998;
t1052 = t944 * t971;
t1051 = (-pkin(3) * t999 - pkin(2) * cos(t1037) - t1026 * pkin(1)) * t969;
t1050 = (-pkin(3) * t1000 - pkin(2) * cos(t1038) - t1027 * pkin(1)) * t970;
t1049 = (-pkin(3) * t1001 - pkin(2) * cos(t1039) - t1028 * pkin(1)) * t971;
t1010 = qJ(1,1) - t1022;
t1009 = qJ(1,1) + t1022;
t1008 = qJ(1,2) - t1021;
t1007 = qJ(1,2) + t1021;
t1006 = qJ(1,3) - t1020;
t1005 = qJ(1,3) + t1020;
t965 = -t1013 * g(1) - t1016 * g(2);
t964 = -t1012 * g(1) - t1015 * g(2);
t963 = -t1011 * g(1) - t1014 * g(2);
t953 = g(3) * t1028 + t968 * t1025;
t952 = g(3) * t1027 + t967 * t1024;
t951 = g(3) * t1026 + t966 * t1023;
t938 = -t962 * pkin(3) + (cos(t995) - cos(t994)) * pkin(2) + (cos(t1010) - cos(t1009)) * pkin(1);
t937 = -t961 * pkin(3) + (cos(t993) - cos(t992)) * pkin(2) + (cos(t1008) - cos(t1007)) * pkin(1);
t936 = -t960 * pkin(3) + (cos(t991) - cos(t990)) * pkin(2) + (cos(t1006) - cos(t1005)) * pkin(1);
t935 = -t959 * pkin(3) + (-sin(t994) - sin(t995)) * pkin(2) + (-sin(t1009) - sin(t1010)) * pkin(1);
t934 = -t958 * pkin(3) + (-sin(t992) - sin(t993)) * pkin(2) + (-sin(t1007) - sin(t1008)) * pkin(1);
t933 = -t957 * pkin(3) + (-sin(t990) - sin(t991)) * pkin(2) + (-sin(t1005) - sin(t1006)) * pkin(1);
t1 = [0, t1064 / 0.2e1 + t1062 / 0.2e1 + t1060 / 0.2e1, t951 * t1066 + t952 * t1067 + t953 * t1068, t1011 * t963 + t1012 * t964 + t1013 * t965 + (t1060 + t1062 + t1064) * t1058, 0, t939 * t1066 + t940 * t1067 + t941 * t1068 + (t935 * t1055 + t934 * t1056 + t933 * t1057) * t1072, t942 * t1066 + t943 * t1067 + t944 * t1068 + (t935 * t1052 + t934 * t1053 + t933 * t1054) * t1072, -g(1); 0, t1065 / 0.2e1 + t1063 / 0.2e1 + t1061 / 0.2e1, t951 * t1069 + t952 * t1070 + t953 * t1071, t1014 * t963 + t1015 * t964 + t1016 * t965 + (t1061 + t1063 + t1065) * t1058, 0, t939 * t1069 + t940 * t1070 + t941 * t1071 + (t938 * t1055 + t937 * t1056 + t936 * t1057) * t1072, t942 * t1069 + t943 * t1070 + t944 * t1071 + (t938 * t1052 + t937 * t1053 + t936 * t1054) * t1072, -g(2); 0, t1059, t953 * t1040 + t952 * t1041 + t951 * t1042, t1059 * pkin(1), 0, t940 * t1041 + t941 * t1040 + t939 * t1042 + (t941 * t1049 + t940 * t1050 + t939 * t1051) * t1030, t943 * t1041 + t944 * t1040 + t942 * t1042 + (t944 * t1049 + t943 * t1050 + t942 * t1051) * t1030, -g(3);];
tau_reg  = t1;
