% Calculate inertia matrix for parallel robot
% P3RRRRR7V2G3A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 10:47
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR7V2G3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:20:32
% EndTime: 2020-08-07 10:20:36
% DurationCPUTime: 4.03s
% Computational Cost: add. (9966->462), mult. (16596->654), div. (792->14), fcn. (8658->84), ass. (0->321)
t1333 = rSges(2,3) + pkin(5);
t1332 = 2 * pkin(1);
t1331 = m(2) / 0.2e1;
t1330 = m(3) / 0.2e1;
t1329 = Icges(2,2) / 0.2e1;
t1328 = Icges(3,2) / 0.2e1;
t1177 = pkin(2) ^ 2;
t1327 = t1177 / 0.2e1;
t1130 = qJ(2,1) + qJ(3,1);
t1086 = sin(t1130);
t1301 = qJ(3,1) + qJ(1,1);
t1094 = qJ(2,1) + t1301;
t1302 = -qJ(3,1) + qJ(1,1);
t1095 = -qJ(2,1) + t1302;
t1169 = 2 * qJ(2,1);
t1129 = t1169 + qJ(3,1);
t1143 = sin(qJ(2,1));
t1144 = sin(qJ(1,1));
t1168 = 2 * qJ(3,1);
t1296 = qJ(1,1) - (2 * qJ(2,1));
t1142 = sin(qJ(3,1));
t1307 = t1142 * pkin(1);
t1318 = pkin(6) + pkin(5);
t1115 = (pkin(7) + t1318);
t1323 = 2 * t1115;
t1207 = ((-sin(t1095) - sin(t1094)) * t1332 + (cos(t1095) + cos(t1094)) * t1323 + (sin((2 * qJ(3,1)) - t1296) - sin(qJ(1,1) + t1169 + t1168) - 0.2e1 * t1144) * pkin(3) + (sin(qJ(3,1) - t1296) - sin(qJ(1,1) + t1129) - sin(t1302) - sin(t1301)) * pkin(2)) / (-t1177 * sin((qJ(2,1) - qJ(3,1))) + pkin(2) * (0.2e1 * t1307 + pkin(2) * t1086 + (sin((t1168 + qJ(2,1))) - t1143) * pkin(3))) / 0.2e1;
t1128 = qJ(2,2) + qJ(3,2);
t1085 = sin(t1128);
t1299 = qJ(3,2) + qJ(1,2);
t1092 = qJ(2,2) + t1299;
t1300 = -qJ(3,2) + qJ(1,2);
t1093 = -qJ(2,2) + t1300;
t1166 = 2 * qJ(2,2);
t1127 = qJ(3,2) + t1166;
t1140 = sin(qJ(2,2));
t1141 = sin(qJ(1,2));
t1165 = 2 * qJ(3,2);
t1295 = qJ(1,2) - (2 * qJ(2,2));
t1139 = sin(qJ(3,2));
t1309 = t1139 * pkin(1);
t1208 = ((-sin(t1093) - sin(t1092)) * t1332 + (cos(t1093) + cos(t1092)) * t1323 + (sin((2 * qJ(3,2)) - t1295) - sin(qJ(1,2) + t1166 + t1165) - 0.2e1 * t1141) * pkin(3) + (sin(qJ(3,2) - t1295) - sin(qJ(1,2) + t1127) - sin(t1300) - sin(t1299)) * pkin(2)) / (-t1177 * sin((qJ(2,2) - qJ(3,2))) + pkin(2) * (0.2e1 * t1309 + pkin(2) * t1085 + (sin((t1165 + qJ(2,2))) - t1140) * pkin(3))) / 0.2e1;
t1126 = qJ(2,3) + qJ(3,3);
t1084 = sin(t1126);
t1297 = qJ(3,3) + qJ(1,3);
t1090 = qJ(2,3) + t1297;
t1298 = -qJ(3,3) + qJ(1,3);
t1091 = -qJ(2,3) + t1298;
t1163 = 2 * qJ(2,3);
t1125 = t1163 + qJ(3,3);
t1137 = sin(qJ(2,3));
t1138 = sin(qJ(1,3));
t1162 = 2 * qJ(3,3);
t1294 = qJ(1,3) - (2 * qJ(2,3));
t1136 = sin(qJ(3,3));
t1311 = t1136 * pkin(1);
t1209 = ((-sin(t1091) - sin(t1090)) * t1332 + (cos(t1091) + cos(t1090)) * t1323 + (sin((2 * qJ(3,3)) - t1294) - sin(qJ(1,3) + t1163 + t1162) - 0.2e1 * t1138) * pkin(3) + (sin(qJ(3,3) - t1294) - sin(qJ(1,3) + t1125) - sin(t1298) - sin(t1297)) * pkin(2)) / (-t1177 * sin((qJ(2,3) - qJ(3,3))) + pkin(2) * (0.2e1 * t1311 + pkin(2) * t1084 + (sin((t1162 + qJ(2,3))) - t1137) * pkin(3))) / 0.2e1;
t1326 = -2 * pkin(1);
t1324 = 0.2e1 * pkin(3);
t1146 = cos(qJ(2,3));
t1322 = 0.2e1 * t1146 ^ 2;
t1149 = cos(qJ(2,2));
t1321 = 0.2e1 * t1149 ^ 2;
t1152 = cos(qJ(2,1));
t1320 = 0.2e1 * t1152 ^ 2;
t1319 = pkin(2) * m(3);
t1317 = m(2) * rSges(2,1);
t1316 = m(2) * rSges(2,2);
t1315 = m(3) * rSges(3,2);
t1172 = rSges(2,2) ^ 2;
t1174 = rSges(2,1) ^ 2;
t1314 = m(3) * t1327 + (-t1172 + t1174) * t1331 + t1329 - Icges(2,1) / 0.2e1;
t1171 = rSges(3,2) ^ 2;
t1173 = rSges(3,1) ^ 2;
t1313 = (-t1171 + t1173) * t1330 - Icges(3,1) / 0.2e1 + t1328;
t1113 = rSges(3,3) + t1318;
t1312 = m(3) * t1113;
t1145 = cos(qJ(3,3));
t1119 = t1145 ^ 2;
t1105 = pkin(3) * t1119;
t1148 = cos(qJ(3,2));
t1121 = t1148 ^ 2;
t1106 = pkin(3) * t1121;
t1151 = cos(qJ(3,1));
t1123 = t1151 ^ 2;
t1107 = pkin(3) * t1123;
t1310 = t1136 * pkin(3);
t1308 = t1139 * pkin(3);
t1306 = t1142 * pkin(3);
t1305 = t1145 * pkin(2);
t1102 = t1145 * pkin(3);
t1304 = t1148 * pkin(2);
t1103 = t1148 * pkin(3);
t1303 = t1151 * pkin(2);
t1104 = t1151 * pkin(3);
t1233 = t1171 + t1173;
t1062 = t1233 * m(3) + Icges(3,3);
t1147 = cos(qJ(1,3));
t1237 = pkin(1) * t1147 + t1138 * t1115;
t1243 = t1136 * t1137;
t1016 = t1237 * t1243 + (t1119 - 0.1e1) * t1147 * pkin(3);
t1133 = legFrame(3,2);
t1096 = sin(t1133);
t1289 = t1016 * t1096;
t1099 = cos(t1133);
t1288 = t1016 * t1099;
t1150 = cos(qJ(1,2));
t1236 = pkin(1) * t1150 + t1141 * t1115;
t1242 = t1139 * t1140;
t1017 = t1236 * t1242 + (t1121 - 0.1e1) * t1150 * pkin(3);
t1134 = legFrame(2,2);
t1097 = sin(t1134);
t1287 = t1017 * t1097;
t1100 = cos(t1134);
t1286 = t1017 * t1100;
t1153 = cos(qJ(1,1));
t1235 = pkin(1) * t1153 + t1144 * t1115;
t1241 = t1142 * t1143;
t1018 = t1235 * t1241 + (t1123 - 0.1e1) * t1153 * pkin(3);
t1135 = legFrame(1,2);
t1098 = sin(t1135);
t1285 = t1018 * t1098;
t1101 = cos(t1135);
t1284 = t1018 * t1101;
t1060 = rSges(3,2) * t1312 - Icges(3,6);
t1061 = rSges(3,1) * t1312 - Icges(3,5);
t1087 = cos(t1126);
t1025 = -t1060 * t1087 - t1084 * t1061;
t1176 = 0.1e1 / pkin(3);
t1283 = t1025 * t1176;
t1088 = cos(t1128);
t1026 = -t1060 * t1088 - t1085 * t1061;
t1282 = t1026 * t1176;
t1089 = cos(t1130);
t1027 = -t1060 * t1089 - t1086 * t1061;
t1281 = t1027 * t1176;
t1239 = rSges(3,1) * t1319;
t1070 = t1145 * t1239;
t1238 = pkin(2) * t1315;
t1191 = t1136 * t1238;
t1034 = t1070 - t1191 + t1062;
t1280 = t1034 * t1176;
t1071 = t1148 * t1239;
t1190 = t1139 * t1238;
t1035 = t1071 - t1190 + t1062;
t1279 = t1035 * t1176;
t1072 = t1151 * t1239;
t1189 = t1142 * t1238;
t1036 = t1072 - t1189 + t1062;
t1278 = t1036 * t1176;
t1225 = pkin(3) * t1243;
t1074 = t1102 + pkin(2);
t1259 = t1074 * t1146;
t1277 = 0.1e1 / (pkin(1) - t1225 + t1259) / t1136;
t1224 = pkin(3) * t1242;
t1075 = t1103 + pkin(2);
t1256 = t1075 * t1149;
t1276 = 0.1e1 / (pkin(1) - t1224 + t1256) / t1139;
t1223 = pkin(3) * t1241;
t1076 = t1104 + pkin(2);
t1253 = t1076 * t1152;
t1275 = 0.1e1 / (pkin(1) - t1223 + t1253) / t1142;
t1045 = 0.1e1 / (t1146 * pkin(2) + pkin(3) * t1087 + pkin(1));
t1274 = t1045 * t1138;
t1273 = t1045 * t1147;
t1046 = 0.1e1 / (t1149 * pkin(2) + pkin(3) * t1088 + pkin(1));
t1272 = t1046 * t1141;
t1271 = t1046 * t1150;
t1047 = 0.1e1 / (t1152 * pkin(2) + pkin(3) * t1089 + pkin(1));
t1270 = t1047 * t1144;
t1269 = t1047 * t1153;
t1055 = pkin(1) * t1137 - t1310;
t1268 = t1055 * t1145;
t1056 = pkin(1) * t1140 - t1308;
t1267 = t1056 * t1148;
t1057 = pkin(1) * t1143 - t1306;
t1266 = t1057 * t1151;
t1265 = t1062 * t1176;
t1161 = pkin(2) / 0.2e1;
t1264 = (t1102 + t1161) * t1136;
t1263 = (t1103 + t1161) * t1139;
t1262 = (t1104 + t1161) * t1142;
t1261 = t1074 * t1096;
t1260 = t1074 * t1099;
t1258 = t1075 * t1097;
t1257 = t1075 * t1100;
t1255 = t1076 * t1098;
t1254 = t1076 * t1101;
t1252 = t1096 * t1147;
t1251 = t1097 * t1150;
t1250 = t1098 * t1153;
t1249 = t1099 * t1147;
t1248 = t1100 * t1150;
t1247 = t1101 * t1153;
t1175 = pkin(3) ^ 2;
t1246 = t1119 * t1175;
t1245 = t1121 * t1175;
t1244 = t1123 * t1175;
t1178 = 0.1e1 / pkin(2);
t1240 = t1176 * t1178;
t1234 = -t1175 / 0.2e1 + t1327;
t1232 = t1172 + t1174;
t1231 = pkin(2) * t1102;
t1230 = pkin(2) * t1103;
t1229 = pkin(2) * t1104;
t1228 = t1074 * t1310;
t1227 = t1075 * t1308;
t1226 = t1076 * t1306;
t1111 = -t1175 + t1177;
t1022 = pkin(1) * t1310 + (t1111 + 0.2e1 * t1231 + 0.2e1 * t1246) * t1137;
t1194 = t1147 * t1243;
t1031 = t1194 * t1324 - t1237;
t1041 = t1231 + t1234 + t1246;
t992 = (t1041 * t1252 - t1099 * t1228) * t1322 + (-t1099 * t1022 - t1031 * t1261) * t1146 - pkin(3) * t1289 - t1055 * t1260;
t1222 = t992 * t1277;
t993 = (-t1041 * t1249 - t1096 * t1228) * t1322 + (-t1096 * t1022 + t1031 * t1260) * t1146 + pkin(3) * t1288 - t1055 * t1261;
t1221 = t993 * t1277;
t1028 = t1311 + (-pkin(3) + t1305 + 0.2e1 * t1105) * t1137;
t1160 = -pkin(3) / 0.2e1;
t1049 = t1105 + t1305 / 0.2e1 + t1160;
t1183 = pkin(2) * t1194 + t1031 * t1145;
t998 = (-t1049 * t1252 + t1099 * t1264) * t1322 + (t1099 * t1028 + t1183 * t1096) * t1146 + t1289 + t1099 * t1268;
t1220 = t998 * t1277;
t999 = (t1049 * t1249 + t1096 * t1264) * t1322 + (t1096 * t1028 - t1183 * t1099) * t1146 - t1288 + t1096 * t1268;
t1219 = t999 * t1277;
t1023 = pkin(1) * t1308 + (t1111 + 0.2e1 * t1230 + 0.2e1 * t1245) * t1140;
t1193 = t1150 * t1242;
t1032 = t1193 * t1324 - t1236;
t1042 = t1230 + t1234 + t1245;
t994 = (t1042 * t1251 - t1100 * t1227) * t1321 + (-t1100 * t1023 - t1032 * t1258) * t1149 - pkin(3) * t1287 - t1056 * t1257;
t1218 = t994 * t1276;
t995 = (-t1042 * t1248 - t1097 * t1227) * t1321 + (-t1097 * t1023 + t1032 * t1257) * t1149 + pkin(3) * t1286 - t1056 * t1258;
t1217 = t995 * t1276;
t1024 = pkin(1) * t1306 + (t1111 + 0.2e1 * t1229 + 0.2e1 * t1244) * t1143;
t1192 = t1153 * t1241;
t1033 = t1192 * t1324 - t1235;
t1043 = t1229 + t1234 + t1244;
t996 = (t1043 * t1250 - t1101 * t1226) * t1320 + (-t1101 * t1024 - t1033 * t1255) * t1152 - pkin(3) * t1285 - t1057 * t1254;
t1216 = t996 * t1275;
t997 = (-t1043 * t1247 - t1098 * t1226) * t1320 + (-t1098 * t1024 + t1033 * t1254) * t1152 + pkin(3) * t1284 - t1057 * t1255;
t1215 = t997 * t1275;
t1214 = t1333 * t1316 - Icges(2,6);
t1029 = t1309 + (-pkin(3) + t1304 + 0.2e1 * t1106) * t1140;
t1050 = t1106 + t1304 / 0.2e1 + t1160;
t1182 = pkin(2) * t1193 + t1032 * t1148;
t1000 = (-t1050 * t1251 + t1100 * t1263) * t1321 + (t1100 * t1029 + t1182 * t1097) * t1149 + t1287 + t1100 * t1267;
t1213 = t1000 * t1276;
t1001 = (t1050 * t1248 + t1097 * t1263) * t1321 + (t1097 * t1029 - t1182 * t1100) * t1149 - t1286 + t1097 * t1267;
t1212 = t1001 * t1276;
t1030 = t1307 + (-pkin(3) + t1303 + 0.2e1 * t1107) * t1143;
t1051 = t1107 + t1303 / 0.2e1 + t1160;
t1181 = pkin(2) * t1192 + t1033 * t1151;
t1002 = (-t1051 * t1250 + t1101 * t1262) * t1320 + (t1101 * t1030 + t1181 * t1098) * t1152 + t1285 + t1101 * t1266;
t1211 = t1002 * t1275;
t1003 = (t1051 * t1247 + t1098 * t1262) * t1320 + (t1098 * t1030 - t1181 * t1101) * t1152 - t1284 + t1098 * t1266;
t1210 = t1003 * t1275;
t1007 = t1041 * t1138 * t1322 + ((pkin(1) - 0.2e1 * t1225) * t1138 - t1147 * t1115) * t1259 - pkin(3) * ((pkin(1) * t1243 - pkin(3) + t1105) * t1138 - t1115 * t1194);
t1206 = t1007 * t1277;
t1008 = t1042 * t1141 * t1321 + ((pkin(1) - 0.2e1 * t1224) * t1141 - t1150 * t1115) * t1256 - pkin(3) * ((pkin(1) * t1242 - pkin(3) + t1106) * t1141 - t1115 * t1193);
t1205 = t1008 * t1276;
t1009 = t1043 * t1144 * t1320 + ((pkin(1) - 0.2e1 * t1223) * t1144 - t1153 * t1115) * t1253 - pkin(3) * ((pkin(1) * t1241 - pkin(3) + t1107) * t1144 - t1115 * t1192);
t1204 = t1009 * t1275;
t1203 = t1178 * t1277;
t1202 = t1178 * t1276;
t1201 = t1178 * t1275;
t1200 = t1096 * t1274;
t1199 = t1099 * t1274;
t1198 = t1097 * t1272;
t1197 = t1100 * t1272;
t1196 = t1098 * t1270;
t1195 = t1101 * t1270;
t1188 = t1177 + t1233;
t1187 = t1232 * m(2) + t1188 * m(3) + Icges(2,3) + Icges(3,3);
t1186 = t1007 * t1176 * t1203;
t1185 = t1008 * t1176 * t1202;
t1184 = t1009 * t1176 * t1201;
t1180 = -pkin(2) * t1312 - t1333 * t1317 + Icges(2,5);
t1159 = 2 * pkin(1) ^ 2;
t1179 = Icges(1,3) + ((2 * pkin(5) ^ 2) + t1159 + ((4 * pkin(5) + 2 * rSges(2,3)) * rSges(2,3)) + t1232) * t1331 + (0.2e1 * t1113 ^ 2 + t1159 + t1188) * t1330 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t1328 + t1329 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t1083 = -rSges(2,1) * t1316 + Icges(2,4);
t1082 = -rSges(3,1) * t1315 + Icges(3,4);
t1081 = 0.2e1 * t1130;
t1080 = 0.2e1 * t1128;
t1079 = 0.2e1 * t1126;
t1073 = t1317 + t1319;
t1021 = 0.2e1 * t1072 + t1187 - 0.2e1 * t1189;
t1020 = 0.2e1 * t1071 + t1187 - 0.2e1 * t1190;
t1019 = 0.2e1 * t1070 + t1187 - 0.2e1 * t1191;
t1012 = (t1060 * t1142 - t1061 * t1151 + t1180) * t1143 - (t1060 * t1151 + t1061 * t1142 + t1214) * t1152;
t1011 = (t1060 * t1139 - t1061 * t1148 + t1180) * t1140 - (t1060 * t1148 + t1061 * t1139 + t1214) * t1149;
t1010 = (t1060 * t1136 - t1061 * t1145 + t1180) * t1137 - (t1060 * t1145 + t1061 * t1136 + t1214) * t1146;
t991 = cos(t1081) * t1313 + cos(t1169) * t1314 + t1072 + (t1073 * t1152 - t1143 * t1316) * t1332 + t1082 * sin(t1081) + t1083 * sin(t1169) + ((cos(t1129) * pkin(2) + t1089 * t1332) * rSges(3,1) + (t1086 * t1326 + (-sin(t1129) - t1142) * pkin(2)) * rSges(3,2)) * m(3) + t1179;
t990 = cos(t1080) * t1313 + cos(t1166) * t1314 + t1071 + (t1073 * t1149 - t1140 * t1316) * t1332 + t1082 * sin(t1080) + t1083 * sin(t1166) + ((cos(t1127) * pkin(2) + t1088 * t1332) * rSges(3,1) + (t1085 * t1326 + (-sin(t1127) - t1139) * pkin(2)) * rSges(3,2)) * m(3) + t1179;
t989 = cos(t1079) * t1313 + cos(t1163) * t1314 + t1070 + (t1073 * t1146 - t1137 * t1316) * t1332 + t1082 * sin(t1079) + t1083 * sin(t1163) + ((cos(t1125) * pkin(2) + t1087 * t1332) * rSges(3,1) + (t1084 * t1326 + (-sin(t1125) - t1136) * pkin(2)) * rSges(3,2)) * m(3) + t1179;
t988 = -t1027 * t1269 + t1036 * t1207 + t1062 * t1184;
t987 = -t1026 * t1271 + t1035 * t1208 + t1062 * t1185;
t986 = -t1025 * t1273 + t1034 * t1209 + t1062 * t1186;
t985 = -t1012 * t1269 + t1021 * t1207 + t1036 * t1184;
t984 = -t1011 * t1271 + t1020 * t1208 + t1035 * t1185;
t983 = -t1010 * t1273 + t1019 * t1209 + t1034 * t1186;
t982 = -t1027 * t1195 + (t1003 * t1036 + t997 * t1265) * t1201;
t981 = -t1026 * t1197 + (t1001 * t1035 + t995 * t1265) * t1202;
t980 = -t1025 * t1199 + (t1034 * t999 + t993 * t1265) * t1203;
t979 = t1027 * t1196 + (t1002 * t1036 + t996 * t1265) * t1201;
t978 = t1026 * t1198 + (t1000 * t1035 + t994 * t1265) * t1202;
t977 = t1025 * t1200 + (t1034 * t998 + t992 * t1265) * t1203;
t976 = t1012 * t1207 + t1027 * t1184 - t991 * t1269;
t975 = t1011 * t1208 + t1026 * t1185 - t990 * t1271;
t974 = t1010 * t1209 + t1025 * t1186 - t989 * t1273;
t973 = -t1012 * t1195 + (t1003 * t1021 + t997 * t1278) * t1201;
t972 = -t1011 * t1197 + (t1001 * t1020 + t995 * t1279) * t1202;
t971 = -t1010 * t1199 + (t1019 * t999 + t993 * t1280) * t1203;
t970 = t1012 * t1196 + (t1002 * t1021 + t996 * t1278) * t1201;
t969 = t1011 * t1198 + (t1000 * t1020 + t994 * t1279) * t1202;
t968 = t1010 * t1200 + (t1019 * t998 + t992 * t1280) * t1203;
t967 = -t991 * t1195 + (t1003 * t1012 + t997 * t1281) * t1201;
t966 = -t990 * t1197 + (t1001 * t1011 + t995 * t1282) * t1202;
t965 = -t989 * t1199 + (t1010 * t999 + t993 * t1283) * t1203;
t964 = t991 * t1196 + (t1002 * t1012 + t996 * t1281) * t1201;
t963 = t990 * t1198 + (t1000 * t1011 + t994 * t1282) * t1202;
t962 = t989 * t1200 + (t1010 * t998 + t992 * t1283) * t1203;
t1 = [-t965 * t1199 - t966 * t1197 - t967 * t1195 + m(4) + (t972 * t1212 + t973 * t1210 + t971 * t1219 + (t982 * t1215 + t981 * t1217 + t980 * t1221) * t1176) * t1178, t965 * t1200 + t966 * t1198 + t967 * t1196 + (t972 * t1213 + t973 * t1211 + t971 * t1220 + (t982 * t1216 + t981 * t1218 + t980 * t1222) * t1176) * t1178, -t965 * t1273 - t966 * t1271 - t967 * t1269 + (t1204 * t982 + t1205 * t981 + t1206 * t980) * t1240 + t973 * t1207 + t972 * t1208 + t971 * t1209; -t962 * t1199 - t963 * t1197 - t964 * t1195 + (t969 * t1212 + t970 * t1210 + t968 * t1219 + (t979 * t1215 + t978 * t1217 + t977 * t1221) * t1176) * t1178, t962 * t1200 + t963 * t1198 + t964 * t1196 + m(4) + (t969 * t1213 + t970 * t1211 + t968 * t1220 + (t979 * t1216 + t978 * t1218 + t977 * t1222) * t1176) * t1178, -t962 * t1273 - t963 * t1271 - t964 * t1269 + (t1204 * t979 + t1205 * t978 + t1206 * t977) * t1240 + t970 * t1207 + t969 * t1208 + t968 * t1209; -t974 * t1199 - t975 * t1197 - t976 * t1195 + (t984 * t1212 + t985 * t1210 + t983 * t1219 + (t988 * t1215 + t987 * t1217 + t986 * t1221) * t1176) * t1178, t974 * t1200 + t975 * t1198 + t976 * t1196 + (t984 * t1213 + t985 * t1211 + t983 * t1220 + (t1216 * t988 + t1218 * t987 + t1222 * t986) * t1176) * t1178, -t974 * t1273 - t975 * t1271 - t976 * t1269 + m(4) + (t1204 * t988 + t1205 * t987 + t1206 * t986) * t1240 + t985 * t1207 + t984 * t1208 + t983 * t1209;];
MX  = t1;
