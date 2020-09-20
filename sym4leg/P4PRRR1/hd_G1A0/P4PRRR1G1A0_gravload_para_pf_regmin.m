% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4PRRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4x11]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4PRRR1G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1A0_gravload_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1A0_gravload_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRR1G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1A0_gravload_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:14:55
% EndTime: 2020-03-02 20:14:57
% DurationCPUTime: 1.62s
% Computational Cost: add. (1677->153), mult. (1313->279), div. (168->6), fcn. (1440->26), ass. (0->140)
t1194 = legFrame(4,3);
t1180 = sin(t1194);
t1184 = cos(t1194);
t1142 = t1180 * g(1) - t1184 * g(2);
t1146 = t1184 * g(1) + t1180 * g(2);
t1190 = pkin(7) + qJ(2,4);
t1170 = qJ(3,4) + t1190;
t1152 = sin(t1170);
t1153 = cos(t1170);
t1090 = t1142 * t1153 + t1146 * t1152;
t1168 = sin(t1190);
t1169 = cos(t1190);
t1244 = 0.1e1 / (t1152 * t1169 - t1168 * t1153);
t1256 = t1090 * t1244;
t1091 = -t1142 * t1152 + t1146 * t1153;
t1255 = t1091 * t1244;
t1195 = legFrame(3,3);
t1181 = sin(t1195);
t1185 = cos(t1195);
t1143 = t1181 * g(1) - t1185 * g(2);
t1147 = t1185 * g(1) + t1181 * g(2);
t1191 = pkin(7) + qJ(2,3);
t1177 = qJ(3,3) + t1191;
t1158 = sin(t1177);
t1161 = cos(t1177);
t1092 = t1143 * t1161 + t1147 * t1158;
t1171 = sin(t1191);
t1174 = cos(t1191);
t1243 = 0.1e1 / (t1158 * t1174 - t1171 * t1161);
t1254 = t1092 * t1243;
t1196 = legFrame(2,3);
t1182 = sin(t1196);
t1186 = cos(t1196);
t1144 = t1182 * g(1) - t1186 * g(2);
t1148 = t1186 * g(1) + t1182 * g(2);
t1192 = pkin(7) + qJ(2,2);
t1178 = qJ(3,2) + t1192;
t1159 = sin(t1178);
t1162 = cos(t1178);
t1093 = t1144 * t1162 + t1148 * t1159;
t1172 = sin(t1192);
t1175 = cos(t1192);
t1242 = 0.1e1 / (t1159 * t1175 - t1172 * t1162);
t1253 = t1093 * t1242;
t1197 = legFrame(1,3);
t1183 = sin(t1197);
t1187 = cos(t1197);
t1145 = t1183 * g(1) - t1187 * g(2);
t1149 = t1187 * g(1) + t1183 * g(2);
t1193 = pkin(7) + qJ(2,1);
t1179 = qJ(3,1) + t1193;
t1160 = sin(t1179);
t1163 = cos(t1179);
t1094 = t1145 * t1163 + t1149 * t1160;
t1173 = sin(t1193);
t1176 = cos(t1193);
t1241 = 0.1e1 / (t1160 * t1176 - t1173 * t1163);
t1252 = t1094 * t1241;
t1095 = -t1143 * t1158 + t1147 * t1161;
t1251 = t1095 * t1243;
t1096 = -t1144 * t1159 + t1148 * t1162;
t1250 = t1096 * t1242;
t1097 = -t1145 * t1160 + t1149 * t1163;
t1249 = t1097 * t1241;
t1248 = koppelP(1,2);
t1247 = koppelP(2,2);
t1246 = koppelP(3,2);
t1245 = koppelP(4,2);
t1198 = xP(4);
t1188 = sin(t1198);
t1189 = cos(t1198);
t1199 = koppelP(4,1);
t1134 = t1188 * t1199 + t1189 * t1245;
t1138 = -t1188 * t1245 + t1189 * t1199;
t1098 = t1134 * t1184 - t1180 * t1138;
t1102 = t1180 * t1134 + t1138 * t1184;
t1212 = t1098 * t1153 - t1102 * t1152;
t1240 = t1212 * t1244;
t1200 = koppelP(3,1);
t1135 = t1188 * t1200 + t1189 * t1246;
t1139 = -t1188 * t1246 + t1189 * t1200;
t1099 = t1135 * t1185 - t1181 * t1139;
t1103 = t1181 * t1135 + t1139 * t1185;
t1211 = t1099 * t1161 - t1103 * t1158;
t1239 = t1211 * t1243;
t1201 = koppelP(2,1);
t1136 = t1188 * t1201 + t1189 * t1247;
t1140 = -t1188 * t1247 + t1189 * t1201;
t1100 = t1136 * t1186 - t1182 * t1140;
t1104 = t1182 * t1136 + t1140 * t1186;
t1210 = t1100 * t1162 - t1104 * t1159;
t1238 = t1210 * t1242;
t1202 = koppelP(1,1);
t1137 = t1188 * t1202 + t1189 * t1248;
t1141 = -t1188 * t1248 + t1189 * t1202;
t1101 = t1137 * t1187 - t1183 * t1141;
t1105 = t1183 * t1137 + t1141 * t1187;
t1209 = t1101 * t1163 - t1105 * t1160;
t1237 = t1209 * t1241;
t1126 = t1152 * t1180 - t1184 * t1153;
t1220 = t1244 * t1126;
t1208 = t1184 * t1152 + t1180 * t1153;
t1219 = t1244 * t1208;
t1128 = t1158 * t1181 - t1185 * t1161;
t1218 = t1243 * t1128;
t1207 = t1185 * t1158 + t1181 * t1161;
t1217 = t1243 * t1207;
t1129 = t1159 * t1182 - t1186 * t1162;
t1216 = t1242 * t1129;
t1206 = t1186 * t1159 + t1182 * t1162;
t1215 = t1242 * t1206;
t1130 = t1160 * t1183 - t1187 * t1163;
t1214 = t1241 * t1130;
t1205 = t1187 * t1160 + t1183 * t1163;
t1213 = t1241 * t1205;
t1204 = 0.1e1 / pkin(2);
t1203 = 0.1e1 / pkin(3);
t1151 = t1189 * g(1) + t1188 * g(2);
t1150 = t1188 * g(1) - t1189 * g(2);
t1113 = -t1145 * t1173 + t1149 * t1176;
t1112 = -t1144 * t1172 + t1148 * t1175;
t1111 = -t1143 * t1171 + t1147 * t1174;
t1110 = t1145 * t1176 + t1149 * t1173;
t1109 = t1144 * t1175 + t1148 * t1172;
t1108 = t1143 * t1174 + t1147 * t1171;
t1107 = -t1142 * t1168 + t1146 * t1169;
t1106 = t1142 * t1169 + t1146 * t1168;
t1089 = -pkin(2) * (t1173 * t1183 - t1187 * t1176) - t1130 * pkin(3);
t1088 = -pkin(2) * (t1172 * t1182 - t1186 * t1175) - t1129 * pkin(3);
t1087 = -pkin(2) * (t1171 * t1181 - t1185 * t1174) - t1128 * pkin(3);
t1086 = pkin(2) * (t1187 * t1173 + t1183 * t1176) + t1205 * pkin(3);
t1085 = pkin(2) * (t1186 * t1172 + t1182 * t1175) + t1206 * pkin(3);
t1084 = pkin(2) * (t1185 * t1171 + t1181 * t1174) + t1207 * pkin(3);
t1083 = -pkin(2) * (t1168 * t1180 - t1184 * t1169) - t1126 * pkin(3);
t1082 = pkin(2) * (t1184 * t1168 + t1180 * t1169) + t1208 * pkin(3);
t1077 = pkin(2) * (t1101 * t1176 - t1105 * t1173) + t1209 * pkin(3);
t1076 = pkin(2) * (t1100 * t1175 - t1104 * t1172) + t1210 * pkin(3);
t1075 = pkin(2) * (t1099 * t1174 - t1103 * t1171) + t1211 * pkin(3);
t1074 = pkin(2) * (t1098 * t1169 - t1102 * t1168) + t1212 * pkin(3);
t1 = [0, 0, (-t1106 * t1220 - t1108 * t1218 - t1109 * t1216 - t1110 * t1214) * t1204, (-t1107 * t1220 - t1111 * t1218 - t1112 * t1216 - t1113 * t1214) * t1204, 0, (-t1090 * t1220 - t1092 * t1218 - t1093 * t1216 - t1094 * t1214 + (-t1083 * t1256 - t1087 * t1254 - t1088 * t1253 - t1089 * t1252) * t1203) * t1204, (-t1091 * t1220 - t1095 * t1218 - t1096 * t1216 - t1097 * t1214 + (-t1083 * t1255 - t1087 * t1251 - t1088 * t1250 - t1089 * t1249) * t1203) * t1204, 0, 0, 0, -t1188 * t1150 - t1189 * t1151; 0, 0, (t1106 * t1219 + t1108 * t1217 + t1109 * t1215 + t1110 * t1213) * t1204, (t1107 * t1219 + t1111 * t1217 + t1112 * t1215 + t1113 * t1213) * t1204, 0, (t1090 * t1219 + t1092 * t1217 + t1093 * t1215 + t1094 * t1213 + (-t1082 * t1256 - t1084 * t1254 - t1085 * t1253 - t1086 * t1252) * t1203) * t1204, (t1091 * t1219 + t1095 * t1217 + t1096 * t1215 + t1097 * t1213 + (-t1082 * t1255 - t1084 * t1251 - t1085 * t1250 - t1086 * t1249) * t1203) * t1204, 0, 0, 0, t1189 * t1150 - t1188 * t1151; -4 * g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, (-t1106 * t1240 - t1108 * t1239 - t1109 * t1238 - t1110 * t1237) * t1204, (-t1107 * t1240 - t1111 * t1239 - t1112 * t1238 - t1113 * t1237) * t1204, 0, (-t1212 * t1256 - t1211 * t1254 - t1210 * t1253 - t1209 * t1252 + (t1074 * t1256 + t1075 * t1254 + t1076 * t1253 + t1077 * t1252) * t1203) * t1204, (-t1212 * t1255 - t1211 * t1251 - t1210 * t1250 - t1209 * t1249 + (t1074 * t1255 + t1075 * t1251 + t1076 * t1250 + t1077 * t1249) * t1203) * t1204, 0, t1150, t1151, 0;];
tau_reg  = t1;
