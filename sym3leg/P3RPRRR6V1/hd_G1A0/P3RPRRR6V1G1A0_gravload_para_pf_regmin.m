% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR6V1G1A0
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
% tau_reg [3x12]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR6V1G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:32:24
% EndTime: 2020-08-06 18:32:25
% DurationCPUTime: 0.75s
% Computational Cost: add. (833->166), mult. (769->201), div. (60->7), fcn. (576->82), ass. (0->133)
t1240 = 2 * pkin(2);
t1180 = sin(qJ(3,3));
t1194 = 0.2e1 * qJ(3,3);
t1107 = 0.1e1 / (pkin(3) * sin(t1194) + t1180 * t1240 + (sin(pkin(7) + qJ(3,3)) + sin(-pkin(7) + qJ(3,3))) * pkin(1));
t1182 = sin(qJ(3,2));
t1195 = 0.2e1 * qJ(3,2);
t1108 = 0.1e1 / (pkin(3) * sin(t1195) + t1182 * t1240 + (sin(pkin(7) + qJ(3,2)) + sin(-pkin(7) + qJ(3,2))) * pkin(1));
t1184 = sin(qJ(3,1));
t1196 = 0.2e1 * qJ(3,1);
t1109 = 0.1e1 / (pkin(3) * sin(t1196) + t1184 * t1240 + (sin(pkin(7) + qJ(3,1)) + sin(-pkin(7) + qJ(3,1))) * pkin(1));
t1176 = qJ(1,1) + pkin(7);
t1179 = legFrame(1,3);
t1151 = t1179 + t1176;
t1146 = qJ(3,1) + t1151;
t1147 = -qJ(3,1) + t1151;
t1218 = sin(t1146) + sin(t1147);
t1175 = qJ(1,2) + pkin(7);
t1178 = legFrame(2,3);
t1150 = t1178 + t1175;
t1142 = qJ(3,2) + t1150;
t1143 = -qJ(3,2) + t1150;
t1219 = sin(t1142) + sin(t1143);
t1174 = qJ(1,3) + pkin(7);
t1177 = legFrame(3,3);
t1149 = t1177 + t1174;
t1138 = qJ(3,3) + t1149;
t1139 = -qJ(3,3) + t1149;
t1220 = sin(t1138) + sin(t1139);
t1243 = -t1220 * t1107 - t1219 * t1108 - t1218 * t1109;
t1242 = -0.2e1 * pkin(1);
t1241 = -2 * pkin(2);
t1193 = (pkin(6) + pkin(5));
t1239 = -2 * t1193;
t1238 = 2 * t1193;
t1237 = g(3) / 0.2e1;
t1236 = cos(pkin(7)) * pkin(1) + pkin(2);
t1166 = qJ(1,1) + t1179;
t1165 = qJ(1,2) + t1178;
t1164 = qJ(1,3) + t1177;
t1186 = cos(qJ(3,3));
t1168 = sin(t1177);
t1171 = cos(t1177);
t1110 = t1168 * g(1) - t1171 * g(2);
t1113 = t1171 * g(1) + t1168 * g(2);
t1158 = sin(t1174);
t1161 = cos(t1174);
t1202 = -t1110 * t1158 + t1113 * t1161;
t1235 = (-g(3) * t1186 + t1180 * t1202) * t1107;
t1234 = (g(3) * t1180 + t1186 * t1202) * t1107;
t1188 = cos(qJ(3,2));
t1169 = sin(t1178);
t1172 = cos(t1178);
t1111 = t1169 * g(1) - t1172 * g(2);
t1114 = t1172 * g(1) + t1169 * g(2);
t1159 = sin(t1175);
t1162 = cos(t1175);
t1201 = -t1111 * t1159 + t1114 * t1162;
t1233 = (-g(3) * t1188 + t1182 * t1201) * t1108;
t1232 = (g(3) * t1182 + t1188 * t1201) * t1108;
t1190 = cos(qJ(3,1));
t1170 = sin(t1179);
t1173 = cos(t1179);
t1112 = t1170 * g(1) - t1173 * g(2);
t1115 = t1173 * g(1) + t1170 * g(2);
t1160 = sin(t1176);
t1163 = cos(t1176);
t1200 = -t1112 * t1160 + t1115 * t1163;
t1231 = (-g(3) * t1190 + t1184 * t1200) * t1109;
t1230 = (g(3) * t1184 + t1190 * t1200) * t1109;
t1116 = 0.1e1 / (t1186 * pkin(3) + t1236);
t1229 = (t1110 * t1161 + t1158 * t1113) * t1116;
t1117 = 0.1e1 / (t1188 * pkin(3) + t1236);
t1228 = (t1111 * t1162 + t1159 * t1114) * t1117;
t1118 = 0.1e1 / (t1190 * pkin(3) + t1236);
t1227 = (t1112 * t1163 + t1160 * t1115) * t1118;
t1131 = sin(t1149);
t1226 = t1116 * t1131;
t1134 = cos(t1149);
t1225 = t1116 * t1134;
t1132 = sin(t1150);
t1224 = t1117 * t1132;
t1135 = cos(t1150);
t1223 = t1117 * t1135;
t1133 = sin(t1151);
t1222 = t1118 * t1133;
t1136 = cos(t1151);
t1221 = t1118 * t1136;
t1217 = cos(t1138) + cos(t1139);
t1216 = cos(t1142) + cos(t1143);
t1215 = cos(t1146) + cos(t1147);
t1214 = 0.2e1 * t1107;
t1213 = 0.2e1 * t1108;
t1212 = 0.2e1 * t1109;
t1211 = t1180 * t1229;
t1210 = t1186 * t1229;
t1209 = t1182 * t1228;
t1208 = t1188 * t1228;
t1207 = t1184 * t1227;
t1206 = t1190 * t1227;
t1181 = sin(qJ(1,3));
t1187 = cos(qJ(1,3));
t1101 = t1110 * t1187 + t1113 * t1181;
t1183 = sin(qJ(1,2));
t1189 = cos(qJ(1,2));
t1102 = t1111 * t1189 + t1114 * t1183;
t1185 = sin(qJ(1,1));
t1191 = cos(qJ(1,1));
t1103 = t1112 * t1191 + t1115 * t1185;
t1199 = -t1101 * t1226 - t1102 * t1224 - t1103 * t1222;
t1198 = t1101 * t1225 + t1102 * t1223 + t1103 * t1221;
t1197 = 0.1e1 / pkin(3);
t1157 = -qJ(3,1) + t1166;
t1156 = qJ(3,1) + t1166;
t1155 = -qJ(3,2) + t1165;
t1154 = qJ(3,2) + t1165;
t1153 = -qJ(3,3) + t1164;
t1152 = qJ(3,3) + t1164;
t1148 = -0.2e1 * qJ(3,1) + t1151;
t1145 = t1196 + t1151;
t1144 = -0.2e1 * qJ(3,2) + t1150;
t1141 = t1195 + t1150;
t1140 = -0.2e1 * qJ(3,3) + t1149;
t1137 = t1194 + t1149;
t1106 = -t1112 * t1185 + t1115 * t1191;
t1105 = -t1111 * t1183 + t1114 * t1189;
t1104 = -t1110 * t1181 + t1113 * t1187;
t1091 = t1133 * t1239 + cos(t1166) * t1242 + t1136 * t1241 - t1215 * pkin(3);
t1090 = t1132 * t1239 + cos(t1165) * t1242 + t1135 * t1241 - t1216 * pkin(3);
t1089 = t1131 * t1239 + cos(t1164) * t1242 + t1134 * t1241 - t1217 * pkin(3);
t1088 = t1136 * t1238 + sin(t1166) * t1242 + t1133 * t1241 - t1218 * pkin(3);
t1087 = t1135 * t1238 + sin(t1165) * t1242 + t1132 * t1241 - t1219 * pkin(3);
t1086 = t1134 * t1238 + sin(t1164) * t1242 + t1131 * t1241 - t1220 * pkin(3);
t1 = [0, t1199, -t1104 * t1226 - t1105 * t1224 - t1106 * t1222, t1199 * pkin(1) + ((-t1215 * t1212 - t1216 * t1213 - t1217 * t1214) * pkin(2) + (-(cos(t1157) + cos(t1156)) * t1212 - (cos(t1155) + cos(t1154)) * t1213 - (cos(t1153) + cos(t1152)) * t1214) * pkin(1) + t1243 * t1238 + (-(cos(t1148) + cos(t1145) + 0.2e1 * t1136) * t1109 - (cos(t1144) + cos(t1141) + 0.2e1 * t1135) * t1108 - (cos(t1140) + cos(t1137) + 0.2e1 * t1134) * t1107) * pkin(3)) * t1237, 0, 0, 0, 0, 0, -t1131 * t1210 - t1132 * t1208 - t1133 * t1206 + (t1089 * t1235 + t1090 * t1233 + t1091 * t1231) * t1197, t1131 * t1211 + t1132 * t1209 + t1133 * t1207 + (t1089 * t1234 + t1090 * t1232 + t1091 * t1230) * t1197, -g(1); 0, t1198, t1104 * t1225 + t1105 * t1223 + t1106 * t1221, t1198 * pkin(1) + ((-(sin(t1157) + sin(t1156)) * t1212 - (sin(t1155) + sin(t1154)) * t1213 - (sin(t1153) + sin(t1152)) * t1214) * pkin(1) + (t1217 * t1107 + t1216 * t1108 + t1215 * t1109) * t1238 + (-(sin(t1148) + sin(t1145) + 0.2e1 * t1133) * t1109 - (sin(t1144) + sin(t1141) + 0.2e1 * t1132) * t1108 - (sin(t1140) + sin(t1137) + 0.2e1 * t1131) * t1107) * pkin(3) + t1243 * t1240) * t1237, 0, 0, 0, 0, 0, t1134 * t1210 + t1135 * t1208 + t1136 * t1206 + (t1086 * t1235 + t1087 * t1233 + t1088 * t1231) * t1197, -t1134 * t1211 - t1135 * t1209 - t1136 * t1207 + (t1086 * t1234 + t1087 * t1232 + t1088 * t1230) * t1197, -g(2); 0, 0, 0, -0.3e1 * g(3), 0, 0, 0, 0, 0, 0, 0, -g(3);];
tau_reg  = t1;
