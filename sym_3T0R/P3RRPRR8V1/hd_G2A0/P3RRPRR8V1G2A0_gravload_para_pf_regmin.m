% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x13]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR8V1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:59:38
% EndTime: 2020-08-06 19:59:39
% DurationCPUTime: 0.68s
% Computational Cost: add. (687->135), mult. (1287->247), div. (117->9), fcn. (1218->23), ass. (0->125)
t1156 = sin(qJ(1,3));
t1162 = cos(qJ(1,3));
t1152 = legFrame(3,2);
t1139 = sin(t1152);
t1142 = cos(t1152);
t1180 = g(1) * t1142 - g(2) * t1139;
t1111 = g(3) * t1156 - t1180 * t1162;
t1158 = sin(qJ(1,2));
t1164 = cos(qJ(1,2));
t1153 = legFrame(2,2);
t1140 = sin(t1153);
t1143 = cos(t1153);
t1179 = g(1) * t1143 - g(2) * t1140;
t1112 = g(3) * t1158 - t1179 * t1164;
t1160 = sin(qJ(1,1));
t1166 = cos(qJ(1,1));
t1154 = legFrame(1,2);
t1141 = sin(t1154);
t1144 = cos(t1154);
t1178 = g(1) * t1144 - g(2) * t1141;
t1113 = g(3) * t1160 - t1178 * t1166;
t1161 = cos(qJ(2,3));
t1229 = pkin(1) * t1161;
t1163 = cos(qJ(2,2));
t1228 = pkin(1) * t1163;
t1165 = cos(qJ(2,1));
t1227 = pkin(1) * t1165;
t1226 = pkin(2) * cos(qJ(2,3) + pkin(5));
t1225 = pkin(2) * cos(qJ(2,2) + pkin(5));
t1224 = pkin(2) * cos(qJ(2,1) + pkin(5));
t1223 = pkin(2) * sin(pkin(5));
t1132 = t1162 * t1229;
t1099 = g(3) * (-qJ(3,3) * t1162 + t1156 * t1229) - t1180 * (qJ(3,3) * t1156 + t1132);
t1135 = cos(pkin(5)) * pkin(2) + pkin(1);
t1155 = sin(qJ(2,3));
t1177 = t1135 * t1161 - t1155 * t1223;
t1117 = 0.1e1 / t1177;
t1219 = t1099 * t1117;
t1133 = t1164 * t1228;
t1100 = g(3) * (-qJ(3,2) * t1164 + t1158 * t1228) - t1179 * (qJ(3,2) * t1158 + t1133);
t1157 = sin(qJ(2,2));
t1176 = t1135 * t1163 - t1157 * t1223;
t1118 = 0.1e1 / t1176;
t1218 = t1100 * t1118;
t1134 = t1166 * t1227;
t1101 = g(3) * (-qJ(3,1) * t1166 + t1160 * t1227) - t1178 * (qJ(3,1) * t1160 + t1134);
t1159 = sin(qJ(2,1));
t1175 = t1135 * t1165 - t1159 * t1223;
t1119 = 0.1e1 / t1175;
t1217 = t1101 * t1119;
t1149 = pkin(4) + qJ(3,3);
t1145 = 0.1e1 / t1149;
t1216 = t1117 * t1145;
t1150 = pkin(4) + qJ(3,2);
t1146 = 0.1e1 / t1150;
t1215 = t1118 * t1146;
t1151 = pkin(4) + qJ(3,1);
t1147 = 0.1e1 / t1151;
t1214 = t1119 * t1147;
t1129 = 0.1e1 / (t1226 + t1229);
t1213 = t1129 * t1139;
t1212 = t1129 * t1142;
t1130 = 0.1e1 / (t1225 + t1228);
t1211 = t1130 * t1140;
t1210 = t1130 * t1143;
t1131 = 0.1e1 / (t1224 + t1227);
t1209 = t1131 * t1141;
t1208 = t1131 * t1144;
t1207 = t1135 * t1156;
t1206 = t1135 * t1158;
t1205 = t1135 * t1160;
t1204 = t1145 * t1162;
t1203 = t1146 * t1164;
t1202 = t1147 * t1166;
t1201 = t1156 * t1223;
t1200 = t1158 * t1223;
t1199 = t1160 * t1223;
t1198 = t1111 * t1145 * t1161;
t1197 = t1112 * t1146 * t1163;
t1196 = t1113 * t1147 * t1165;
t1195 = t1111 * t1216;
t1194 = t1112 * t1215;
t1193 = t1113 * t1214;
t1192 = t1113 * t1202;
t1174 = g(3) * t1162 + t1180 * t1156;
t1191 = t1174 * t1216;
t1173 = g(3) * t1164 + t1179 * t1158;
t1190 = t1173 * t1215;
t1172 = g(3) * t1166 + t1178 * t1160;
t1189 = t1172 * t1214;
t1188 = t1111 * t1204;
t1187 = t1112 * t1203;
t1186 = t1117 * t1198;
t1185 = t1118 * t1197;
t1184 = t1119 * t1196;
t1183 = t1155 * t1195;
t1182 = t1157 * t1194;
t1181 = t1159 * t1193;
t1126 = g(1) * t1139 + g(2) * t1142;
t1093 = -t1126 * t1161 + t1155 * t1174;
t1127 = g(1) * t1140 + g(2) * t1143;
t1095 = -t1127 * t1163 + t1157 * t1173;
t1128 = g(1) * t1141 + g(2) * t1144;
t1097 = -t1128 * t1165 + t1159 * t1172;
t1171 = t1093 * t1213 + t1095 * t1211 + t1097 * t1209;
t1170 = t1093 * t1212 + t1095 * t1210 + t1097 * t1208;
t1169 = t1172 * t1202 + t1173 * t1203 + t1174 * t1204;
t1087 = (-t1139 * t1207 + t1142 * t1223) * t1161 + t1155 * (t1135 * t1142 + t1139 * t1201);
t1088 = (-t1140 * t1206 + t1143 * t1223) * t1163 + t1157 * (t1135 * t1143 + t1140 * t1200);
t1089 = (-t1141 * t1205 + t1144 * t1223) * t1165 + t1159 * (t1135 * t1144 + t1141 * t1199);
t1168 = t1087 * t1191 + t1088 * t1190 + t1089 * t1189;
t1090 = (t1139 * t1223 + t1142 * t1207) * t1161 + (t1135 * t1139 - t1142 * t1201) * t1155;
t1091 = (t1140 * t1223 + t1143 * t1206) * t1163 + (t1135 * t1140 - t1143 * t1200) * t1157;
t1092 = (t1141 * t1223 + t1144 * t1205) * t1165 + (t1135 * t1141 - t1144 * t1199) * t1159;
t1167 = t1090 * t1191 + t1091 * t1190 + t1092 * t1189;
t1122 = t1135 * t1159 + t1165 * t1223;
t1121 = t1135 * t1157 + t1163 * t1223;
t1120 = t1135 * t1155 + t1161 * t1223;
t1104 = -t1151 * t1166 + t1175 * t1160;
t1103 = -t1150 * t1164 + t1176 * t1158;
t1102 = -t1149 * t1162 + t1177 * t1156;
t1098 = t1128 * t1159 + t1165 * t1172;
t1096 = t1127 * t1157 + t1163 * t1173;
t1094 = t1126 * t1155 + t1161 * t1174;
t1 = [0, t1090 * t1195 + t1091 * t1194 + t1092 * t1193, t1167, 0, 0, 0, 0, 0, t1090 * t1186 + t1091 * t1185 + t1092 * t1184 + t1171, -t1090 * t1183 - t1091 * t1182 - t1092 * t1181 + t1094 * t1213 + t1096 * t1211 + t1098 * t1209, -t1167, (t1092 * t1217 - (t1104 * t1144 + t1122 * t1141) * t1113) * t1147 + (t1091 * t1218 - (t1103 * t1143 + t1121 * t1140) * t1112) * t1146 + (t1090 * t1219 - (t1102 * t1142 + t1120 * t1139) * t1111) * t1145 + t1171 * pkin(1), -g(1); 0, t1087 * t1195 + t1088 * t1194 + t1089 * t1193, t1168, 0, 0, 0, 0, 0, t1087 * t1186 + t1088 * t1185 + t1089 * t1184 + t1170, -t1087 * t1183 - t1088 * t1182 - t1089 * t1181 + t1094 * t1212 + t1096 * t1210 + t1098 * t1208, -t1168, (t1089 * t1217 - (-t1104 * t1141 + t1122 * t1144) * t1113) * t1147 + (t1088 * t1218 - (-t1103 * t1140 + t1121 * t1143) * t1112) * t1146 + (t1087 * t1219 - (-t1102 * t1139 + t1120 * t1142) * t1111) * t1145 + t1170 * pkin(1), -g(2); 0, t1187 + t1188 + t1192, t1169, 0, 0, 0, 0, 0, t1162 * t1198 + t1164 * t1197 + t1166 * t1196, -t1155 * t1188 - t1157 * t1187 - t1159 * t1192, -t1169, (t1166 * t1101 - (t1160 * t1151 + t1166 * t1224 + t1134) * t1113) * t1147 + (t1164 * t1100 - (t1150 * t1158 + t1164 * t1225 + t1133) * t1112) * t1146 + (t1162 * t1099 - (t1149 * t1156 + t1162 * t1226 + t1132) * t1111) * t1145, -g(3);];
tau_reg  = t1;
