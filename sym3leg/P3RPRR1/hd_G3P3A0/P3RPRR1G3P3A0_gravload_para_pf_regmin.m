% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRR1G3P3A0
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
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRR1G3P3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:27:04
% EndTime: 2020-03-09 21:27:05
% DurationCPUTime: 0.66s
% Computational Cost: add. (933->139), mult. (651->170), div. (81->4), fcn. (570->63), ass. (0->114)
t1032 = sin(qJ(1,1));
t1035 = cos(qJ(1,1));
t1029 = legFrame(1,2);
t1021 = sin(t1029);
t1024 = cos(t1029);
t976 = t1024 * g(1) - t1021 * g(2);
t1038 = g(3) * t1035 + t976 * t1032;
t1031 = sin(qJ(1,2));
t1034 = cos(qJ(1,2));
t1028 = legFrame(2,2);
t1020 = sin(t1028);
t1023 = cos(t1028);
t975 = t1023 * g(1) - t1020 * g(2);
t1039 = g(3) * t1034 + t975 * t1031;
t1042 = pkin(7) + qJ(3,1);
t1012 = qJ(1,1) + t1042;
t1006 = sin(t1012);
t979 = 0.1e1 / (pkin(1) * sin(t1042) + sin(qJ(3,1)) * pkin(2));
t1047 = t1006 * t979;
t1041 = pkin(7) + qJ(3,2);
t1011 = qJ(1,2) + t1041;
t1005 = sin(t1011);
t978 = 0.1e1 / (pkin(1) * sin(t1041) + sin(qJ(3,2)) * pkin(2));
t1048 = t1005 * t978;
t1040 = pkin(7) + qJ(3,3);
t1010 = qJ(1,3) + t1040;
t1004 = sin(t1010);
t977 = 0.1e1 / (pkin(1) * sin(t1040) + sin(qJ(3,3)) * pkin(2));
t1049 = t1004 * t977;
t1030 = sin(qJ(1,3));
t1033 = cos(qJ(1,3));
t1027 = legFrame(3,2);
t1019 = sin(t1027);
t1022 = cos(t1027);
t974 = t1022 * g(1) - t1019 * g(2);
t956 = g(3) * t1033 + t974 * t1030;
t1082 = -t1038 * t1047 - t1039 * t1048 - t956 * t1049;
t1037 = 0.1e1 / pkin(3);
t1081 = t1037 / 0.2e1;
t1045 = qJ(1,1) + pkin(7);
t1002 = t1029 + t1045;
t996 = qJ(3,1) + t1002;
t1003 = -t1029 + t1045;
t997 = qJ(3,1) + t1003;
t970 = cos(t997) + cos(t996);
t1053 = t970 * t979;
t1080 = t1053 / 0.2e1;
t1044 = qJ(1,2) + pkin(7);
t1000 = t1028 + t1044;
t994 = qJ(3,2) + t1000;
t1001 = -t1028 + t1044;
t995 = qJ(3,2) + t1001;
t969 = cos(t995) + cos(t994);
t1054 = t969 * t978;
t1079 = t1054 / 0.2e1;
t1043 = qJ(1,3) + pkin(7);
t998 = t1027 + t1043;
t992 = qJ(3,3) + t998;
t999 = -t1027 + t1043;
t993 = qJ(3,3) + t999;
t968 = cos(t993) + cos(t992);
t1055 = t968 * t977;
t1078 = t1055 / 0.2e1;
t1050 = sin(t996) - sin(t997);
t1056 = t1050 * t979;
t1077 = -t1056 / 0.2e1;
t1051 = sin(t994) - sin(t995);
t1057 = t1051 * t978;
t1076 = -t1057 / 0.2e1;
t1052 = sin(t992) - sin(t993);
t1058 = t1052 * t977;
t1075 = -t1058 / 0.2e1;
t1074 = t1038 * t1053;
t1073 = t1038 * t1056;
t1072 = t1039 * t1054;
t1071 = t1039 * t1057;
t1069 = pkin(1) / 0.2e1;
t1007 = cos(t1010);
t947 = g(3) * t1007 + t974 * t1004;
t1067 = t947 * t977;
t948 = -g(3) * t1004 + t974 * t1007;
t1066 = t948 * t977;
t1008 = cos(t1011);
t949 = g(3) * t1008 + t975 * t1005;
t1065 = t949 * t978;
t950 = -g(3) * t1005 + t975 * t1008;
t1064 = t950 * t978;
t1009 = cos(t1012);
t951 = g(3) * t1009 + t976 * t1006;
t1063 = t951 * t979;
t952 = -g(3) * t1006 + t976 * t1009;
t1062 = t952 * t979;
t1061 = (pkin(3) * t1004 + pkin(2) * sin(t1043) + t1030 * pkin(1)) * t977;
t1060 = (pkin(3) * t1005 + pkin(2) * sin(t1044) + t1031 * pkin(1)) * t978;
t1059 = (pkin(3) * t1006 + pkin(2) * sin(t1045) + t1032 * pkin(1)) * t979;
t1018 = qJ(1,1) - t1029;
t1017 = qJ(1,1) + t1029;
t1016 = qJ(1,2) - t1028;
t1015 = qJ(1,2) + t1028;
t1014 = qJ(1,3) - t1027;
t1013 = qJ(1,3) + t1027;
t973 = -t1021 * g(1) - t1024 * g(2);
t972 = -t1020 * g(1) - t1023 * g(2);
t971 = -t1019 * g(1) - t1022 * g(2);
t961 = -g(3) * t1032 + t976 * t1035;
t959 = -g(3) * t1031 + t975 * t1034;
t957 = -g(3) * t1030 + t974 * t1033;
t946 = -t970 * pkin(3) + (-cos(t1003) - cos(t1002)) * pkin(2) + (-cos(t1018) - cos(t1017)) * pkin(1);
t945 = -t969 * pkin(3) + (-cos(t1001) - cos(t1000)) * pkin(2) + (-cos(t1016) - cos(t1015)) * pkin(1);
t944 = -t968 * pkin(3) + (-cos(t999) - cos(t998)) * pkin(2) + (-cos(t1014) - cos(t1013)) * pkin(1);
t943 = t1050 * pkin(3) + (sin(t1002) - sin(t1003)) * pkin(2) + (sin(t1017) - sin(t1018)) * pkin(1);
t942 = t1051 * pkin(3) + (sin(t1000) - sin(t1001)) * pkin(2) + (sin(t1015) - sin(t1016)) * pkin(1);
t941 = t1052 * pkin(3) + (sin(t998) - sin(t999)) * pkin(2) + (sin(t1013) - sin(t1014)) * pkin(1);
t1 = [0, t1074 / 0.2e1 + t1072 / 0.2e1 + t956 * t1078, t957 * t1078 + t959 * t1079 + t961 * t1080, t1019 * t971 + t1020 * t972 + t1021 * t973 + (t956 * t1055 + t1072 + t1074) * t1069, 0, t947 * t1078 + t949 * t1079 + t951 * t1080 + (t946 * t1063 + t945 * t1065 + t944 * t1067) * t1081, t948 * t1078 + t950 * t1079 + t952 * t1080 + (t946 * t1062 + t945 * t1064 + t944 * t1066) * t1081, -g(1); 0, -t1073 / 0.2e1 - t1071 / 0.2e1 + t956 * t1075, t957 * t1075 + t959 * t1076 + t961 * t1077, t1022 * t971 + t1023 * t972 + t1024 * t973 + (-t956 * t1058 - t1071 - t1073) * t1069, 0, t947 * t1075 + t949 * t1076 + t951 * t1077 + (t943 * t1063 + t942 * t1065 + t941 * t1067) * t1081, t948 * t1075 + t950 * t1076 + t952 * t1077 + (t943 * t1062 + t942 * t1064 + t941 * t1066) * t1081, -g(2); 0, t1082, -t961 * t1047 - t1048 * t959 - t1049 * t957, t1082 * pkin(1), 0, -t947 * t1049 - t949 * t1048 - t951 * t1047 + (t951 * t1059 + t949 * t1060 + t947 * t1061) * t1037, -t948 * t1049 - t950 * t1048 - t952 * t1047 + (t952 * t1059 + t950 * t1060 + t948 * t1061) * t1037, -g(3);];
tau_reg  = t1;
