% Calculate Gravitation load for parallel robot
% P3RRPRR8V2G1A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-07 13:12
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:11:03
% EndTime: 2022-11-07 13:11:04
% DurationCPUTime: 0.64s
% Computational Cost: add. (666->126), mult. (1071->186), div. (27->6), fcn. (615->35), ass. (0->95)
t1120 = -m(3) / 0.2e1;
t1119 = pkin(5) + qJ(3,2);
t1118 = pkin(5) + qJ(3,3);
t1117 = qJ(3,1) + pkin(5);
t1065 = qJ(2,3) + pkin(7);
t1116 = pkin(3) * cos(t1065);
t1066 = qJ(2,2) + pkin(7);
t1115 = pkin(3) * cos(t1066);
t1067 = qJ(2,1) + pkin(7);
t1114 = pkin(3) * cos(t1067);
t1068 = sin(pkin(7));
t1112 = cos(pkin(7));
t1011 = -m(2) * rSges(2,1) + (-rSges(3,1) * t1112 + rSges(3,2) * t1068 - pkin(2)) * m(3);
t1113 = g(3) * t1011;
t1111 = 0.2e1 * pkin(2) * pkin(3);
t1110 = 2 * pkin(1);
t1069 = legFrame(3,3);
t1047 = sin(t1069);
t1050 = cos(t1069);
t1026 = -t1047 * g(1) + t1050 * g(2);
t1029 = t1050 * g(1) + t1047 * g(2);
t1062 = -pkin(6) - t1118;
t1053 = 0.1e1 / t1062;
t1073 = sin(qJ(1,3));
t1079 = cos(qJ(1,3));
t1109 = (t1026 * t1079 - t1029 * t1073) * t1053;
t1070 = legFrame(2,3);
t1048 = sin(t1070);
t1051 = cos(t1070);
t1027 = -t1048 * g(1) + t1051 * g(2);
t1030 = t1051 * g(1) + t1048 * g(2);
t1063 = -pkin(6) - t1119;
t1054 = 0.1e1 / t1063;
t1075 = sin(qJ(1,2));
t1081 = cos(qJ(1,2));
t1108 = (t1027 * t1081 - t1030 * t1075) * t1054;
t1071 = legFrame(1,3);
t1049 = sin(t1071);
t1052 = cos(t1071);
t1028 = -t1049 * g(1) + t1052 * g(2);
t1031 = t1052 * g(1) + t1049 * g(2);
t1064 = -pkin(6) - t1117;
t1055 = 0.1e1 / t1064;
t1077 = sin(qJ(1,1));
t1083 = cos(qJ(1,1));
t1107 = (t1028 * t1083 - t1031 * t1077) * t1055;
t1025 = m(2) * rSges(2,2) + (rSges(3,1) * t1068 + t1112 * rSges(3,2)) * m(3);
t1072 = sin(qJ(2,3));
t1078 = cos(qJ(2,3));
t1103 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t1093 = -t1011 * t1078 - t1025 * t1072 + t1103;
t1035 = (-rSges(2,3) - pkin(5)) * m(2) + m(1) * rSges(1,2);
t1102 = -(rSges(3,3) + t1118) * m(3) + t1035;
t1106 = t1053 * ((t1093 * t1073 + t1102 * t1079) * t1029 + (t1102 * t1073 - t1093 * t1079) * t1026);
t1074 = sin(qJ(2,2));
t1080 = cos(qJ(2,2));
t1092 = -t1011 * t1080 - t1025 * t1074 + t1103;
t1101 = -(rSges(3,3) + t1119) * m(3) + t1035;
t1105 = t1054 * ((t1092 * t1075 + t1101 * t1081) * t1030 + (t1101 * t1075 - t1092 * t1081) * t1027);
t1076 = sin(qJ(2,1));
t1082 = cos(qJ(2,1));
t1091 = -t1011 * t1082 - t1025 * t1076 + t1103;
t1100 = -(rSges(3,3) + t1117) * m(3) + t1035;
t1104 = t1055 * ((t1091 * t1077 + t1100 * t1083) * t1031 + (t1100 * t1077 - t1091 * t1083) * t1028);
t1099 = t1072 * pkin(2) + pkin(3) * sin(t1065);
t1098 = t1074 * pkin(2) + pkin(3) * sin(t1066);
t1097 = t1076 * pkin(2) + pkin(3) * sin(t1067);
t1096 = t1026 * t1073 + t1029 * t1079;
t1095 = t1027 * t1075 + t1030 * t1081;
t1094 = t1028 * t1077 + t1031 * t1083;
t1090 = pkin(2) ^ 2;
t1089 = pkin(3) ^ 2;
t1088 = 0.2e1 * qJ(2,1);
t1087 = 0.2e1 * qJ(2,2);
t1086 = 0.2e1 * qJ(2,3);
t1058 = t1082 * pkin(2);
t1057 = t1080 * pkin(2);
t1056 = t1078 * pkin(2);
t1039 = t1058 + pkin(1);
t1038 = t1057 + pkin(1);
t1037 = t1056 + pkin(1);
t1024 = g(3) * t1025;
t1023 = -t1049 * t1077 + t1052 * t1083;
t1022 = t1049 * t1083 + t1052 * t1077;
t1021 = -t1048 * t1075 + t1051 * t1081;
t1020 = t1048 * t1081 + t1051 * t1075;
t1019 = -t1047 * t1073 + t1050 * t1079;
t1018 = t1047 * t1079 + t1050 * t1073;
t1017 = t1039 * t1083 - t1077 * t1064;
t1016 = t1038 * t1081 - t1075 * t1063;
t1015 = t1037 * t1079 - t1073 * t1062;
t1014 = t1077 * t1039 + t1083 * t1064;
t1013 = t1075 * t1038 + t1081 * t1063;
t1012 = t1073 * t1037 + t1079 * t1062;
t1 = [-t1019 * t1106 - t1021 * t1105 - t1023 * t1104 - m(4) * g(1) + (-(-t1049 * t1014 + t1017 * t1052 + t1023 * t1114) * t1107 - (-t1048 * t1013 + t1016 * t1051 + t1021 * t1115) * t1108 - (-t1047 * t1012 + t1015 * t1050 + t1019 * t1116) * t1109) * m(3); -t1018 * t1106 - t1020 * t1105 - t1022 * t1104 - m(4) * g(2) + (-(t1014 * t1052 + t1017 * t1049 + t1022 * t1114) * t1107 - (t1013 * t1051 + t1016 * t1048 + t1020 * t1115) * t1108 - (t1012 * t1050 + t1015 * t1047 + t1018 * t1116) * t1109) * m(3); -m(4) * g(3) + (-t1097 * t1104 + (-t1094 * t1011 + t1024) * t1076 + t1082 * (t1094 * t1025 + t1113) + (sin(t1088 + pkin(7)) * t1111 + t1089 * sin(0.2e1 * t1067) + t1090 * sin(t1088) + t1097 * t1110) * t1107 * t1120) / (t1058 + t1114) + (-t1098 * t1105 + (-t1095 * t1011 + t1024) * t1074 + t1080 * (t1095 * t1025 + t1113) + (sin(t1087 + pkin(7)) * t1111 + t1089 * sin(0.2e1 * t1066) + t1090 * sin(t1087) + t1098 * t1110) * t1108 * t1120) / (t1057 + t1115) + (-t1099 * t1106 + (-t1096 * t1011 + t1024) * t1072 + t1078 * (t1096 * t1025 + t1113) + (sin(t1086 + pkin(7)) * t1111 + t1089 * sin(0.2e1 * t1065) + t1090 * sin(t1086) + t1099 * t1110) * t1109 * t1120) / (t1056 + t1116);];
taugX  = t1;
