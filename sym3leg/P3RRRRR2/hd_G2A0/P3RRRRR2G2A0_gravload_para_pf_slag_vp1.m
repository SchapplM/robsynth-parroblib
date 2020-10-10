% Calculate Gravitation load for parallel robot
% P3RRRRR2G2A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 21:10
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR2G2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G2A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G2A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G2A0_gravload_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR2G2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G2A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:08:41
% EndTime: 2020-03-09 21:08:42
% DurationCPUTime: 1.09s
% Computational Cost: add. (657->161), mult. (1263->256), div. (69->14), fcn. (657->48), ass. (0->134)
t1215 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t1145 = cos(qJ(3,1));
t1214 = t1145 ^ 2;
t1142 = cos(qJ(3,2));
t1213 = t1142 ^ 2;
t1139 = cos(qJ(3,3));
t1212 = t1139 ^ 2;
t1136 = sin(qJ(3,1));
t1211 = t1145 * rSges(3,1) - t1136 * rSges(3,2);
t1133 = sin(qJ(3,2));
t1210 = t1142 * rSges(3,1) - t1133 * rSges(3,2);
t1130 = sin(qJ(3,3));
t1209 = t1139 * rSges(3,1) - t1130 * rSges(3,2);
t1208 = -2 * pkin(1);
t1206 = m(1) * rSges(1,2);
t1205 = rSges(3,3) * g(3);
t1132 = sin(qJ(1,3));
t1204 = pkin(1) * t1132;
t1135 = sin(qJ(1,2));
t1203 = pkin(1) * t1135;
t1138 = sin(qJ(1,1));
t1202 = pkin(1) * t1138;
t1127 = legFrame(3,2);
t1104 = sin(t1127);
t1107 = cos(t1127);
t1078 = t1107 * g(1) - t1104 * g(2);
t1201 = rSges(3,1) * t1078;
t1128 = legFrame(2,2);
t1105 = sin(t1128);
t1108 = cos(t1128);
t1079 = t1108 * g(1) - t1105 * g(2);
t1200 = rSges(3,1) * t1079;
t1129 = legFrame(1,2);
t1106 = sin(t1129);
t1109 = cos(t1129);
t1080 = t1109 * g(1) - t1106 * g(2);
t1199 = rSges(3,1) * t1080;
t1198 = rSges(3,3) * t1078;
t1197 = rSges(3,3) * t1079;
t1196 = rSges(3,3) * t1080;
t1189 = qJ(2,1) - qJ(3,1);
t1188 = qJ(2,1) + qJ(3,1);
t1187 = qJ(2,2) - qJ(3,2);
t1186 = qJ(2,2) + qJ(3,2);
t1185 = qJ(2,3) - qJ(3,3);
t1184 = qJ(2,3) + qJ(3,3);
t1149 = rSges(2,2) * g(3);
t1066 = m(2) * (-rSges(2,1) * t1078 + t1149);
t1151 = rSges(2,1) * g(3);
t1069 = m(2) * (rSges(2,2) * t1078 + t1151);
t1072 = rSges(3,2) * t1078;
t1098 = qJ(1,3) + t1184;
t1085 = cos(t1098);
t1099 = qJ(1,3) + t1185;
t1086 = cos(t1099);
t1124 = qJ(1,3) + qJ(2,3);
t1092 = sin(t1124);
t1095 = cos(t1124);
t1110 = g(3) * t1206;
t1131 = sin(qJ(2,3));
t1112 = 0.1e1 / t1131;
t1141 = cos(qJ(1,3));
t1148 = rSges(3,2) * g(3);
t1150 = rSges(3,1) * g(3);
t1167 = t1215 * g(3);
t1168 = m(3) * t1205;
t1183 = ((t1078 * t1206 + t1167) * t1132 + (-(t1148 + t1201) * t1086 / 0.2e1 - (t1072 - t1150) * sin(t1099) / 0.2e1 + (t1148 - t1201) * t1085 / 0.2e1 + (t1072 + t1150) * sin(t1098) / 0.2e1) * m(3) + (t1066 - t1168) * t1095 + t1092 * (-m(3) * t1198 + t1069) + (-t1078 * t1215 + t1110) * t1141) * t1112;
t1067 = m(2) * (-rSges(2,1) * t1079 + t1149);
t1070 = m(2) * (rSges(2,2) * t1079 + t1151);
t1073 = rSges(3,2) * t1079;
t1100 = qJ(1,2) + t1186;
t1087 = cos(t1100);
t1101 = qJ(1,2) + t1187;
t1088 = cos(t1101);
t1125 = qJ(1,2) + qJ(2,2);
t1093 = sin(t1125);
t1096 = cos(t1125);
t1134 = sin(qJ(2,2));
t1113 = 0.1e1 / t1134;
t1144 = cos(qJ(1,2));
t1182 = ((t1079 * t1206 + t1167) * t1135 + (-(t1148 + t1200) * t1088 / 0.2e1 - (t1073 - t1150) * sin(t1101) / 0.2e1 + (t1148 - t1200) * t1087 / 0.2e1 + (t1073 + t1150) * sin(t1100) / 0.2e1) * m(3) + (t1067 - t1168) * t1096 + t1093 * (-m(3) * t1197 + t1070) + (-t1079 * t1215 + t1110) * t1144) * t1113;
t1068 = m(2) * (-rSges(2,1) * t1080 + t1149);
t1071 = m(2) * (rSges(2,2) * t1080 + t1151);
t1074 = rSges(3,2) * t1080;
t1102 = qJ(1,1) + t1188;
t1089 = cos(t1102);
t1103 = qJ(1,1) + t1189;
t1090 = cos(t1103);
t1126 = qJ(1,1) + qJ(2,1);
t1094 = sin(t1126);
t1097 = cos(t1126);
t1137 = sin(qJ(2,1));
t1114 = 0.1e1 / t1137;
t1147 = cos(qJ(1,1));
t1181 = ((t1080 * t1206 + t1167) * t1138 + (-(t1148 + t1199) * t1090 / 0.2e1 - (t1074 - t1150) * sin(t1103) / 0.2e1 + (t1148 - t1199) * t1089 / 0.2e1 + (t1074 + t1150) * sin(t1102) / 0.2e1) * m(3) + (t1068 - t1168) * t1097 + t1094 * (-m(3) * t1196 + t1071) + (-t1080 * t1215 + t1110) * t1147) * t1114;
t1116 = 0.1e1 / t1139;
t1180 = (-(t1104 * g(1) + t1107 * g(2)) * t1209 + (g(3) * t1095 + t1092 * t1078) * (t1130 * rSges(3,1) + t1139 * rSges(3,2))) * t1116;
t1119 = 0.1e1 / t1142;
t1179 = (-(t1105 * g(1) + t1108 * g(2)) * t1210 + (g(3) * t1096 + t1093 * t1079) * (t1133 * rSges(3,1) + t1142 * rSges(3,2))) * t1119;
t1122 = 0.1e1 / t1145;
t1178 = (-(t1106 * g(1) + t1109 * g(2)) * t1211 + (g(3) * t1097 + t1094 * t1080) * (t1136 * rSges(3,1) + t1145 * rSges(3,2))) * t1122;
t1140 = cos(qJ(2,3));
t1075 = t1141 * t1131 + t1132 * t1140;
t1177 = t1075 * t1139;
t1143 = cos(qJ(2,2));
t1076 = t1144 * t1134 + t1135 * t1143;
t1176 = t1076 * t1142;
t1146 = cos(qJ(2,1));
t1077 = t1147 * t1137 + t1138 * t1146;
t1175 = t1077 * t1145;
t1174 = t1104 * t1130;
t1173 = t1105 * t1133;
t1172 = t1106 * t1136;
t1171 = t1107 * t1130;
t1170 = t1108 * t1133;
t1169 = t1109 * t1136;
t1166 = t1075 * t1212 * pkin(2);
t1165 = t1076 * t1213 * pkin(2);
t1164 = t1077 * t1214 * pkin(2);
t1163 = t1140 * t1130 * pkin(1);
t1162 = t1143 * t1133 * pkin(1);
t1161 = t1146 * t1136 * pkin(1);
t1159 = t1116 * t1183;
t1158 = t1119 * t1182;
t1157 = t1122 * t1181;
t1060 = t1066 * t1095 + t1092 * t1069 + ((-t1209 * t1078 - t1205) * t1095 + t1092 * (t1209 * g(3) - t1198)) * m(3);
t1156 = t1060 * t1112 / t1212;
t1061 = t1067 * t1096 + t1093 * t1070 + ((-t1210 * t1079 - t1205) * t1096 + t1093 * (t1210 * g(3) - t1197)) * m(3);
t1155 = t1061 * t1113 / t1213;
t1062 = t1068 * t1097 + t1094 * t1071 + ((-t1211 * t1080 - t1205) * t1097 + t1094 * (t1211 * g(3) - t1196)) * m(3);
t1154 = t1062 * t1114 / t1214;
t1153 = 1 / pkin(1);
t1152 = 0.1e1 / pkin(2);
t1 = [-m(4) * g(1) + ((t1109 * t1175 + t1172) * t1157 + (t1108 * t1176 + t1173) * t1158 + (t1107 * t1177 + t1174) * t1159) * t1153 + (((-t1109 * t1164 + (-pkin(2) * t1172 - t1109 * t1202) * t1145 - t1106 * t1161) * t1154 + (-t1108 * t1165 + (-pkin(2) * t1173 - t1108 * t1203) * t1142 - t1105 * t1162) * t1155 + (-t1107 * t1166 + (-pkin(2) * t1174 - t1107 * t1204) * t1139 - t1104 * t1163) * t1156) * t1153 + (t1104 * t1180 + t1105 * t1179 + t1106 * t1178) * m(3)) * t1152; -m(4) * g(2) + ((-t1106 * t1175 + t1169) * t1157 + (-t1105 * t1176 + t1170) * t1158 + (-t1104 * t1177 + t1171) * t1159) * t1153 + (((t1106 * t1164 + (-pkin(2) * t1169 + t1106 * t1202) * t1145 - t1109 * t1161) * t1154 + (t1105 * t1165 + (-pkin(2) * t1170 + t1105 * t1203) * t1142 - t1108 * t1162) * t1155 + (t1104 * t1166 + (-pkin(2) * t1171 + t1104 * t1204) * t1139 - t1107 * t1163) * t1156) * t1153 + (t1107 * t1180 + t1108 * t1179 + t1109 * t1178) * m(3)) * t1152; -m(4) * g(3) + (t1095 * t1183 + t1096 * t1182 + t1097 * t1181 + ((t1147 * t1208 + (-t1089 - t1090) * pkin(2)) / (sin(t1188) + sin(t1189)) * t1062 + (t1144 * t1208 + (-t1087 - t1088) * pkin(2)) / (sin(t1186) + sin(t1187)) * t1061 + (t1141 * t1208 + (-t1085 - t1086) * pkin(2)) / (sin(t1184) + sin(t1185)) * t1060) * t1152) * t1153;];
taugX  = t1;
