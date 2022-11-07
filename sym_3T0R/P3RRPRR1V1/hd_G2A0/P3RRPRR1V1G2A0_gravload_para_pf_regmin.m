% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x15]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-04 17:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR1V1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_regmin: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:08:20
% EndTime: 2022-11-04 17:08:21
% DurationCPUTime: 0.59s
% Computational Cost: add. (588->108), mult. (996->207), div. (177->7), fcn. (1077->18), ass. (0->119)
t1213 = sin(qJ(1,3));
t1219 = cos(qJ(1,3));
t1209 = legFrame(3,2);
t1193 = sin(t1209);
t1196 = cos(t1209);
t1238 = t1196 * g(1) - t1193 * g(2);
t1175 = g(3) * t1213 - t1238 * t1219;
t1215 = sin(qJ(1,2));
t1221 = cos(qJ(1,2));
t1210 = legFrame(2,2);
t1194 = sin(t1210);
t1197 = cos(t1210);
t1237 = t1197 * g(1) - t1194 * g(2);
t1176 = g(3) * t1215 - t1237 * t1221;
t1217 = sin(qJ(1,1));
t1223 = cos(qJ(1,1));
t1211 = legFrame(1,2);
t1195 = sin(t1211);
t1198 = cos(t1211);
t1236 = t1198 * g(1) - t1195 * g(2);
t1177 = g(3) * t1217 - t1236 * t1223;
t1218 = cos(qJ(2,3));
t1259 = t1218 * t1219;
t1264 = t1213 * t1218;
t1166 = g(3) * (pkin(1) * t1264 - t1219 * qJ(3,3)) - t1238 * (pkin(1) * t1259 + t1213 * qJ(3,3));
t1202 = 0.1e1 / t1218;
t1283 = t1166 * t1202;
t1220 = cos(qJ(2,2));
t1257 = t1220 * t1221;
t1262 = t1215 * t1220;
t1167 = g(3) * (pkin(1) * t1262 - t1221 * qJ(3,2)) - t1237 * (pkin(1) * t1257 + t1215 * qJ(3,2));
t1203 = 0.1e1 / t1220;
t1282 = t1167 * t1203;
t1222 = cos(qJ(2,1));
t1255 = t1222 * t1223;
t1260 = t1217 * t1222;
t1168 = g(3) * (pkin(1) * t1260 - t1223 * qJ(3,1)) - t1236 * (pkin(1) * t1255 + t1217 * qJ(3,1));
t1204 = 0.1e1 / t1222;
t1281 = t1168 * t1204;
t1206 = pkin(3) + qJ(3,3);
t1199 = 0.1e1 / t1206;
t1280 = t1175 * t1199;
t1207 = pkin(3) + qJ(3,2);
t1200 = 0.1e1 / t1207;
t1279 = t1176 * t1200;
t1208 = pkin(3) + qJ(3,1);
t1201 = 0.1e1 / t1208;
t1278 = t1177 * t1201;
t1277 = t1193 * t1202;
t1276 = t1194 * t1203;
t1275 = t1195 * t1204;
t1274 = t1196 * t1202;
t1273 = t1197 * t1203;
t1272 = t1198 * t1204;
t1271 = t1199 * t1202;
t1270 = t1199 * t1219;
t1269 = t1200 * t1203;
t1268 = t1200 * t1221;
t1267 = t1201 * t1204;
t1266 = t1201 * t1223;
t1212 = sin(qJ(2,3));
t1224 = pkin(1) + pkin(2);
t1265 = t1212 * t1224;
t1214 = sin(qJ(2,2));
t1263 = t1214 * t1224;
t1216 = sin(qJ(2,1));
t1261 = t1216 * t1224;
t1258 = t1218 * t1224;
t1256 = t1220 * t1224;
t1254 = t1222 * t1224;
t1181 = -t1193 * t1264 + t1212 * t1196;
t1253 = t1181 * t1271;
t1182 = -t1194 * t1262 + t1214 * t1197;
t1252 = t1182 * t1269;
t1183 = -t1195 * t1260 + t1216 * t1198;
t1251 = t1183 * t1267;
t1184 = t1193 * t1212 + t1196 * t1264;
t1250 = t1184 * t1271;
t1185 = t1194 * t1214 + t1197 * t1262;
t1249 = t1185 * t1269;
t1186 = t1195 * t1216 + t1198 * t1260;
t1248 = t1186 * t1267;
t1247 = t1175 * t1270;
t1246 = t1176 * t1268;
t1245 = t1177 * t1266;
t1244 = t1175 * t1253;
t1243 = t1176 * t1252;
t1242 = t1177 * t1251;
t1241 = t1175 * t1250;
t1240 = t1176 * t1249;
t1239 = t1177 * t1248;
t1235 = -t1206 * t1219 + t1213 * t1258;
t1234 = -t1207 * t1221 + t1215 * t1256;
t1233 = -t1208 * t1223 + t1217 * t1254;
t1232 = g(3) * t1219 + t1238 * t1213;
t1231 = g(3) * t1221 + t1237 * t1215;
t1230 = g(3) * t1223 + t1236 * t1217;
t1229 = t1230 * t1266 + t1231 * t1268 + t1232 * t1270;
t1190 = t1193 * g(1) + t1196 * g(2);
t1160 = -t1190 * t1218 + t1212 * t1232;
t1191 = t1194 * g(1) + t1197 * g(2);
t1162 = -t1191 * t1220 + t1214 * t1231;
t1192 = t1195 * g(1) + t1198 * g(2);
t1164 = -t1192 * t1222 + t1216 * t1230;
t1205 = 0.1e1 / t1224;
t1228 = (t1160 * t1277 + t1162 * t1276 + t1164 * t1275) * t1205;
t1227 = (t1160 * t1274 + t1162 * t1273 + t1164 * t1272) * t1205;
t1226 = t1230 * t1251 + t1231 * t1252 + t1232 * t1253;
t1225 = t1230 * t1248 + t1231 * t1249 + t1232 * t1250;
t1165 = t1192 * t1216 + t1222 * t1230;
t1163 = t1191 * t1214 + t1220 * t1231;
t1161 = t1190 * t1212 + t1218 * t1232;
t1159 = -t1212 * t1247 - t1214 * t1246 - t1216 * t1245;
t1158 = t1255 * t1278 + t1257 * t1279 + t1259 * t1280;
t1157 = t1184 * t1280 + t1185 * t1279 + t1186 * t1278 + t1228;
t1156 = t1181 * t1280 + t1182 * t1279 + t1183 * t1278 + t1227;
t1155 = -t1212 * t1241 - t1214 * t1240 - t1216 * t1239 + (t1161 * t1277 + t1163 * t1276 + t1165 * t1275) * t1205;
t1154 = -t1212 * t1244 - t1214 * t1243 - t1216 * t1242 + (t1161 * t1274 + t1163 * t1273 + t1165 * t1272) * t1205;
t1 = [0, t1239 + t1240 + t1241, t1225, 0, 0, 0, 0, 0, t1157, t1155, t1157, t1155, -t1225, (t1186 * t1281 - (t1195 * t1261 + t1233 * t1198) * t1177) * t1201 + (t1185 * t1282 - (t1194 * t1263 + t1234 * t1197) * t1176) * t1200 + (t1184 * t1283 - (t1193 * t1265 + t1235 * t1196) * t1175) * t1199 + pkin(1) * t1228, -g(1); 0, t1242 + t1243 + t1244, t1226, 0, 0, 0, 0, 0, t1156, t1154, t1156, t1154, -t1226, (t1183 * t1281 - (-t1233 * t1195 + t1198 * t1261) * t1177) * t1201 + (t1182 * t1282 - (-t1234 * t1194 + t1197 * t1263) * t1176) * t1200 + (t1181 * t1283 - (-t1235 * t1193 + t1196 * t1265) * t1175) * t1199 + pkin(1) * t1227, -g(2); 0, t1245 + t1246 + t1247, t1229, 0, 0, 0, 0, 0, t1158, t1159, t1158, t1159, -t1229, (t1223 * t1168 - (t1217 * t1208 + t1223 * t1254) * t1177) * t1201 + (t1221 * t1167 - (t1215 * t1207 + t1221 * t1256) * t1176) * t1200 + (t1219 * t1166 - (t1213 * t1206 + t1219 * t1258) * t1175) * t1199, -g(3);];
tau_reg  = t1;
