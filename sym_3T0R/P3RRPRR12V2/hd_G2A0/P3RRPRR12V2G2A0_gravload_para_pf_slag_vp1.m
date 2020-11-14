% Calculate Gravitation load for parallel robot
% P3RRPRR12V2G2A0
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
% Datum: 2020-08-06 19:23
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V2G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:20:41
% EndTime: 2020-08-06 19:20:42
% DurationCPUTime: 1.00s
% Computational Cost: add. (870->181), mult. (1173->317), div. (45->6), fcn. (648->18), ass. (0->133)
t1109 = m(2) * rSges(2,2);
t1065 = (-qJ(3,1) - rSges(3,3)) * m(3) + t1109;
t1066 = (pkin(2) + rSges(3,1)) * m(3) + m(2) * rSges(2,1);
t1101 = sin(qJ(2,1));
t1107 = cos(qJ(2,1));
t1181 = -t1065 * t1101 + t1066 * t1107;
t1064 = (-qJ(3,2) - rSges(3,3)) * m(3) + t1109;
t1099 = sin(qJ(2,2));
t1105 = cos(qJ(2,2));
t1180 = -t1064 * t1099 + t1066 * t1105;
t1063 = (-qJ(3,3) - rSges(3,3)) * m(3) + t1109;
t1097 = sin(qJ(2,3));
t1103 = cos(qJ(2,3));
t1179 = -t1063 * t1097 + t1066 * t1103;
t1111 = (pkin(2) + pkin(3));
t1178 = -2 * t1111;
t1094 = legFrame(3,2);
t1076 = sin(t1094);
t1177 = t1076 * qJ(3,3);
t1095 = legFrame(2,2);
t1077 = sin(t1095);
t1176 = t1077 * qJ(3,2);
t1096 = legFrame(1,2);
t1078 = sin(t1096);
t1175 = t1078 * qJ(3,1);
t1079 = cos(t1094);
t1174 = t1079 * qJ(3,3);
t1080 = cos(t1095);
t1173 = t1080 * qJ(3,2);
t1081 = cos(t1096);
t1172 = t1081 * qJ(3,1);
t1098 = sin(qJ(1,3));
t1171 = t1098 * qJ(3,3);
t1100 = sin(qJ(1,2));
t1170 = t1100 * qJ(3,2);
t1102 = sin(qJ(1,1));
t1169 = t1102 * qJ(3,1);
t1054 = t1076 * g(1) + t1079 * g(2);
t1057 = t1079 * g(1) - t1076 * g(2);
t1104 = cos(qJ(1,3));
t1117 = g(3) * t1104 + t1057 * t1098;
t1032 = (-t1066 * t1054 + t1117 * t1063) * t1103 + (t1054 * t1063 + t1117 * t1066) * t1097;
t1112 = 0.1e1 / qJ(3,3);
t1168 = t1032 * t1112;
t1055 = t1077 * g(1) + t1080 * g(2);
t1058 = t1080 * g(1) - t1077 * g(2);
t1106 = cos(qJ(1,2));
t1116 = g(3) * t1106 + t1058 * t1100;
t1033 = (-t1066 * t1055 + t1116 * t1064) * t1105 + (t1055 * t1064 + t1116 * t1066) * t1099;
t1113 = 0.1e1 / qJ(3,2);
t1167 = t1033 * t1113;
t1056 = t1078 * g(1) + t1081 * g(2);
t1059 = t1081 * g(1) - t1078 * g(2);
t1108 = cos(qJ(1,1));
t1115 = g(3) * t1108 + t1059 * t1102;
t1034 = (-t1066 * t1056 + t1115 * t1065) * t1107 + (t1056 * t1065 + t1115 * t1066) * t1101;
t1114 = 0.1e1 / qJ(3,1);
t1166 = t1034 * t1114;
t1110 = pkin(5) - pkin(6);
t1123 = pkin(1) * t1098 - t1110 * t1104;
t1132 = t1097 * t1171;
t1165 = (t1123 + 0.2e1 * t1132) * t1111;
t1122 = pkin(1) * t1100 - t1110 * t1106;
t1131 = t1099 * t1170;
t1164 = (t1122 + 0.2e1 * t1131) * t1111;
t1121 = pkin(1) * t1102 - t1110 * t1108;
t1130 = t1101 * t1169;
t1163 = (t1121 + 0.2e1 * t1130) * t1111;
t1073 = t1097 * qJ(3,3);
t1138 = t1111 * t1103;
t1120 = t1073 + pkin(1) + t1138;
t1051 = 0.1e1 / t1120;
t1162 = t1051 * t1112;
t1074 = t1099 * qJ(3,2);
t1137 = t1111 * t1105;
t1119 = t1074 + pkin(1) + t1137;
t1052 = 0.1e1 / t1119;
t1161 = t1052 * t1113;
t1075 = t1101 * qJ(3,1);
t1136 = t1111 * t1107;
t1118 = t1075 + pkin(1) + t1136;
t1053 = 0.1e1 / t1118;
t1160 = t1053 * t1114;
t1153 = (qJ(3,3) + t1111) * (-qJ(3,3) + t1111);
t1152 = (qJ(3,2) + t1111) * (-qJ(3,2) + t1111);
t1151 = (qJ(3,1) + t1111) * (-qJ(3,1) + t1111);
t1150 = t1097 * t1111;
t1070 = t1098 * t1110;
t1149 = t1098 * t1111;
t1148 = t1099 * t1111;
t1071 = t1100 * t1110;
t1147 = t1100 * t1111;
t1146 = t1101 * t1111;
t1072 = t1102 * t1110;
t1145 = t1102 * t1111;
t1060 = (-rSges(3,2) - pkin(5)) * m(3) + (-rSges(2,3) - pkin(5)) * m(2) + m(1) * rSges(1,2);
t1050 = t1060 * g(3);
t1062 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t1061 = t1062 * g(3);
t1029 = t1050 * t1104 + (t1179 * g(3) + t1061) * t1098 + ((-t1062 - t1179) * t1104 + t1060 * t1098) * t1057;
t1144 = t1104 * t1029;
t1030 = t1050 * t1106 + (t1180 * g(3) + t1061) * t1100 + ((-t1062 - t1180) * t1106 + t1060 * t1100) * t1058;
t1143 = t1106 * t1030;
t1031 = t1050 * t1108 + (t1181 * g(3) + t1061) * t1102 + ((-t1062 - t1181) * t1108 + t1060 * t1102) * t1059;
t1142 = t1108 * t1031;
t1067 = pkin(1) * t1097 + qJ(3,3);
t1141 = t1111 * t1067;
t1068 = pkin(1) * t1099 + qJ(3,2);
t1140 = t1111 * t1068;
t1069 = pkin(1) * t1101 + qJ(3,1);
t1139 = t1111 * t1069;
t1135 = qJ(3,1) * t1178;
t1134 = qJ(3,2) * t1178;
t1133 = qJ(3,3) * t1178;
t1129 = (-t1054 * t1103 + t1117 * t1097) * t1162;
t1128 = (-t1055 * t1105 + t1116 * t1099) * t1161;
t1127 = (-t1056 * t1107 + t1115 * t1101) * t1160;
t1126 = t1098 * t1153;
t1125 = t1100 * t1152;
t1124 = t1102 * t1151;
t1093 = t1107 ^ 2;
t1092 = t1105 ^ 2;
t1091 = t1103 ^ 2;
t1049 = pkin(1) * qJ(3,1) - t1101 * t1151;
t1048 = pkin(1) * qJ(3,2) - t1099 * t1152;
t1047 = pkin(1) * qJ(3,3) - t1097 * t1153;
t1046 = t1121 + t1130;
t1044 = t1122 + t1131;
t1042 = t1123 + t1132;
t1040 = t1121 * t1101 + t1169;
t1039 = t1122 * t1099 + t1170;
t1038 = t1123 * t1097 + t1171;
t1 = [-m(4) * g(1) + (t1081 * t1142 + ((t1081 * t1145 - t1175) * t1093 + (t1046 * t1081 + t1078 * t1146) * t1107 + t1078 * t1069) * t1166) * t1053 + (t1080 * t1143 + ((t1080 * t1147 - t1176) * t1092 + (t1044 * t1080 + t1077 * t1148) * t1105 + t1077 * t1068) * t1167) * t1052 + (t1079 * t1144 + ((t1079 * t1149 - t1177) * t1091 + (t1042 * t1079 + t1076 * t1150) * t1103 + t1076 * t1067) * t1168) * t1051 + (-((t1078 * t1135 + t1081 * t1124) * t1093 + (-t1078 * t1049 + t1081 * t1163) * t1107 + t1040 * t1172 + t1078 * t1139) * t1127 - ((t1077 * t1134 + t1080 * t1125) * t1092 + (-t1077 * t1048 + t1080 * t1164) * t1105 + t1039 * t1173 + t1077 * t1140) * t1128 - ((t1076 * t1133 + t1079 * t1126) * t1091 + (-t1076 * t1047 + t1079 * t1165) * t1103 + t1038 * t1174 + t1076 * t1141) * t1129) * m(3); -m(4) * g(2) + (-t1078 * t1142 + ((-t1078 * t1145 - t1172) * t1093 + (-t1046 * t1078 + t1081 * t1146) * t1107 + t1081 * t1069) * t1166) * t1053 + (-t1077 * t1143 + ((-t1077 * t1147 - t1173) * t1092 + (-t1044 * t1077 + t1080 * t1148) * t1105 + t1080 * t1068) * t1167) * t1052 + (-t1076 * t1144 + ((-t1076 * t1149 - t1174) * t1091 + (-t1042 * t1076 + t1079 * t1150) * t1103 + t1079 * t1067) * t1168) * t1051 + (-((-t1078 * t1124 + t1081 * t1135) * t1093 + (-t1081 * t1049 - t1078 * t1163) * t1107 - t1040 * t1175 + t1081 * t1139) * t1127 - ((-t1077 * t1125 + t1080 * t1134) * t1092 + (-t1080 * t1048 - t1077 * t1164) * t1105 - t1039 * t1176 + t1080 * t1140) * t1128 - ((-t1076 * t1126 + t1079 * t1133) * t1091 + (-t1079 * t1047 - t1076 * t1165) * t1103 - t1038 * t1177 + t1079 * t1141) * t1129) * m(3); -t1102 * t1053 * t1031 + (t1118 * t1108 + t1072) * t1107 * t1034 * t1160 - t1100 * t1052 * t1030 + (t1119 * t1106 + t1071) * t1105 * t1033 * t1161 - t1098 * t1051 * t1029 + (t1120 * t1104 + t1070) * t1103 * t1032 * t1162 - m(4) * g(3) + (-(t1108 * t1093 * t1151 + ((0.2e1 * t1075 + pkin(1)) * t1108 + t1072) * t1136 + qJ(3,1) * (t1069 * t1108 + t1101 * t1072)) * t1127 - (t1106 * t1092 * t1152 + ((0.2e1 * t1074 + pkin(1)) * t1106 + t1071) * t1137 + qJ(3,2) * (t1068 * t1106 + t1099 * t1071)) * t1128 - (t1104 * t1091 * t1153 + ((0.2e1 * t1073 + pkin(1)) * t1104 + t1070) * t1138 + qJ(3,3) * (t1067 * t1104 + t1097 * t1070)) * t1129) * m(3);];
taugX  = t1;
