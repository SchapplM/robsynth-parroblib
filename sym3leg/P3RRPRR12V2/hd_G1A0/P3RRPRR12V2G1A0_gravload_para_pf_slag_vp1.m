% Calculate Gravitation load for parallel robot
% P3RRPRR12V2G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2020-08-06 19:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V2G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:14:50
% EndTime: 2020-08-06 19:14:50
% DurationCPUTime: 0.57s
% Computational Cost: add. (717->157), mult. (1011->252), div. (36->6), fcn. (636->18), ass. (0->113)
t1100 = pkin(2) + pkin(3);
t1052 = (pkin(2) + rSges(3,1)) * m(3) + m(2) * rSges(2,1);
t1140 = t1052 * g(3);
t1086 = sin(qJ(2,3));
t1139 = t1086 * qJ(3,3);
t1088 = sin(qJ(2,2));
t1138 = t1088 * qJ(3,2);
t1090 = sin(qJ(2,1));
t1137 = t1090 * qJ(3,1);
t1092 = cos(qJ(2,3));
t1101 = 0.1e1 / qJ(3,3);
t1083 = legFrame(3,3);
t1068 = sin(t1083);
t1071 = cos(t1083);
t1038 = -t1068 * g(1) + t1071 * g(2);
t1041 = t1071 * g(1) + t1068 * g(2);
t1087 = sin(qJ(1,3));
t1093 = cos(qJ(1,3));
t1112 = t1038 * t1087 + t1041 * t1093;
t1136 = (-g(3) * t1092 + t1112 * t1086) * t1101;
t1094 = cos(qJ(2,2));
t1102 = 0.1e1 / qJ(3,2);
t1084 = legFrame(2,3);
t1069 = sin(t1084);
t1072 = cos(t1084);
t1039 = -t1069 * g(1) + t1072 * g(2);
t1042 = t1072 * g(1) + t1069 * g(2);
t1089 = sin(qJ(1,2));
t1095 = cos(qJ(1,2));
t1111 = t1039 * t1089 + t1042 * t1095;
t1135 = (-g(3) * t1094 + t1111 * t1088) * t1102;
t1096 = cos(qJ(2,1));
t1103 = 0.1e1 / qJ(3,1);
t1085 = legFrame(1,3);
t1070 = sin(t1085);
t1073 = cos(t1085);
t1040 = -t1070 * g(1) + t1073 * g(2);
t1043 = t1073 * g(1) + t1070 * g(2);
t1091 = sin(qJ(1,1));
t1097 = cos(qJ(1,1));
t1110 = t1040 * t1091 + t1043 * t1097;
t1134 = (-g(3) * t1096 + t1110 * t1090) * t1103;
t1099 = pkin(5) - pkin(6);
t1062 = t1087 * t1099;
t1063 = t1089 * t1099;
t1064 = t1091 * t1099;
t1065 = t1099 * t1093;
t1066 = t1099 * t1095;
t1067 = t1099 * t1097;
t1133 = t1100 * t1092;
t1132 = t1100 * t1094;
t1131 = t1100 * t1096;
t1098 = m(2) * rSges(2,2);
t1049 = (-qJ(3,3) - rSges(3,3)) * m(3) + t1098;
t1130 = t1101 * ((t1112 * t1049 - t1140) * t1092 + (g(3) * t1049 + t1112 * t1052) * t1086);
t1050 = (-qJ(3,2) - rSges(3,3)) * m(3) + t1098;
t1129 = t1102 * ((t1111 * t1050 - t1140) * t1094 + (g(3) * t1050 + t1111 * t1052) * t1088);
t1051 = (-qJ(3,1) - rSges(3,3)) * m(3) + t1098;
t1128 = t1103 * ((t1110 * t1051 - t1140) * t1096 + (g(3) * t1051 + t1110 * t1052) * t1090);
t1127 = t1092 * t1130;
t1126 = t1094 * t1129;
t1125 = t1096 * t1128;
t1054 = pkin(1) + t1139;
t1044 = 0.1e1 / (t1054 + t1133);
t1124 = t1044 * t1136;
t1056 = pkin(1) + t1138;
t1045 = 0.1e1 / (t1056 + t1132);
t1123 = t1045 * t1135;
t1058 = pkin(1) + t1137;
t1046 = 0.1e1 / (t1058 + t1131);
t1122 = t1046 * t1134;
t1121 = (qJ(3,3) + t1100) * (-qJ(3,3) + t1100) * t1092 ^ 2;
t1120 = (qJ(3,2) + t1100) * (-qJ(3,2) + t1100) * t1094 ^ 2;
t1119 = (qJ(3,1) + t1100) * (-qJ(3,1) + t1100) * t1096 ^ 2;
t1053 = pkin(1) + 0.2e1 * t1139;
t1118 = t1053 * t1093 + t1062;
t1117 = t1054 * t1093 + t1062;
t1055 = pkin(1) + 0.2e1 * t1138;
t1116 = t1055 * t1095 + t1063;
t1115 = t1056 * t1095 + t1063;
t1057 = pkin(1) + 0.2e1 * t1137;
t1114 = t1057 * t1097 + t1064;
t1113 = t1058 * t1097 + t1064;
t1059 = pkin(1) * t1086 + qJ(3,3);
t1109 = t1059 * t1093 + t1086 * t1062;
t1060 = pkin(1) * t1088 + qJ(3,2);
t1108 = t1060 * t1095 + t1088 * t1063;
t1061 = pkin(1) * t1090 + qJ(3,1);
t1107 = t1061 * t1097 + t1090 * t1064;
t1048 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t1106 = t1049 * t1086 - t1052 * t1092 - t1048;
t1105 = t1050 * t1088 - t1052 * t1094 - t1048;
t1104 = t1051 * t1090 - t1052 * t1096 - t1048;
t1047 = (-rSges(3,2) - pkin(5)) * m(3) + (-rSges(2,3) - pkin(5)) * m(2) + m(1) * rSges(1,2);
t1037 = t1091 * t1058 - t1067;
t1036 = t1091 * t1057 - t1067;
t1035 = t1089 * t1056 - t1066;
t1034 = t1089 * t1055 - t1066;
t1033 = t1087 * t1054 - t1065;
t1032 = t1087 * t1053 - t1065;
t1031 = t1070 * t1097 + t1073 * t1091;
t1030 = t1070 * t1091 - t1073 * t1097;
t1029 = t1069 * t1095 + t1072 * t1089;
t1028 = t1069 * t1089 - t1072 * t1095;
t1027 = t1068 * t1093 + t1071 * t1087;
t1026 = t1068 * t1087 - t1071 * t1093;
t1025 = t1091 * t1061 - t1090 * t1067;
t1024 = t1089 * t1060 - t1088 * t1066;
t1023 = t1087 * t1059 - t1086 * t1065;
t1016 = (t1047 * t1097 - t1091 * t1104) * t1043 + (t1091 * t1047 + t1104 * t1097) * t1040;
t1015 = (t1047 * t1095 - t1089 * t1105) * t1042 + (t1089 * t1047 + t1105 * t1095) * t1039;
t1014 = (t1047 * t1093 - t1087 * t1106) * t1041 + (t1087 * t1047 + t1106 * t1093) * t1038;
t1 = [-m(4) * g(1) + (-t1031 * t1016 - (t1030 * t1131 + t1037 * t1070 - t1113 * t1073) * t1125) * t1046 + (-t1029 * t1015 - (t1028 * t1132 + t1035 * t1069 - t1115 * t1072) * t1126) * t1045 + (-t1027 * t1014 - (t1026 * t1133 + t1033 * t1068 - t1117 * t1071) * t1127) * t1044 + (-(-t1030 * t1119 - (t1070 * t1036 - t1114 * t1073) * t1131 - (t1070 * t1025 - t1107 * t1073) * qJ(3,1)) * t1122 - (-t1028 * t1120 - (t1069 * t1034 - t1116 * t1072) * t1132 - (t1069 * t1024 - t1108 * t1072) * qJ(3,2)) * t1123 - (-t1026 * t1121 - (t1068 * t1032 - t1118 * t1071) * t1133 - (t1068 * t1023 - t1109 * t1071) * qJ(3,3)) * t1124) * m(3); -m(4) * g(2) + (-t1030 * t1016 + (t1031 * t1131 + t1037 * t1073 + t1070 * t1113) * t1125) * t1046 + (-t1028 * t1015 + (t1029 * t1132 + t1035 * t1072 + t1069 * t1115) * t1126) * t1045 + (-t1026 * t1014 + (t1027 * t1133 + t1033 * t1071 + t1068 * t1117) * t1127) * t1044 + (-(t1031 * t1119 + (t1036 * t1073 + t1114 * t1070) * t1131 + (t1025 * t1073 + t1107 * t1070) * qJ(3,1)) * t1122 - (t1029 * t1120 + (t1034 * t1072 + t1116 * t1069) * t1132 + (t1024 * t1072 + t1108 * t1069) * qJ(3,2)) * t1123 - (t1027 * t1121 + (t1032 * t1071 + t1118 * t1068) * t1133 + (t1023 * t1071 + t1109 * t1068) * qJ(3,3)) * t1124) * m(3); t1086 * t1130 + t1088 * t1129 + t1090 * t1128 - m(4) * g(3) + (-(-t1096 * qJ(3,1) + t1100 * t1090) * t1134 - (-t1094 * qJ(3,2) + t1100 * t1088) * t1135 - (-t1092 * qJ(3,3) + t1100 * t1086) * t1136) * m(3);];
taugX  = t1;
