% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRR1G1P1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:22
% EndTime: 2020-03-09 21:23:23
% DurationCPUTime: 1.22s
% Computational Cost: add. (10696->160), mult. (7576->274), div. (1086->10), fcn. (4824->32), ass. (0->143)
t1157 = pkin(1) ^ 2;
t1224 = MDP(4) * t1157 + MDP(1);
t1139 = pkin(7) + qJ(3,3);
t1147 = sin(qJ(3,3));
t1101 = pkin(1) * sin(t1139) + t1147 * pkin(2);
t1092 = 0.1e1 / t1101;
t1140 = pkin(7) + qJ(3,2);
t1148 = sin(qJ(3,2));
t1102 = pkin(1) * sin(t1140) + t1148 * pkin(2);
t1095 = 0.1e1 / t1102;
t1141 = pkin(7) + qJ(3,1);
t1149 = sin(qJ(3,1));
t1103 = pkin(1) * sin(t1141) + t1149 * pkin(2);
t1098 = 0.1e1 / t1103;
t1093 = 0.1e1 / t1101 ^ 2;
t1096 = 0.1e1 / t1102 ^ 2;
t1099 = 0.1e1 / t1103 ^ 2;
t1142 = sin(pkin(7));
t1223 = pkin(1) * t1142;
t1143 = cos(pkin(7));
t1222 = pkin(2) * t1143;
t1150 = cos(qJ(3,3));
t1221 = pkin(2) * t1150;
t1151 = cos(qJ(3,2));
t1220 = pkin(2) * t1151;
t1152 = cos(qJ(3,1));
t1219 = pkin(2) * t1152;
t1218 = 0.2e1 * pkin(2) * pkin(3);
t1135 = legFrame(3,3) + qJ(1,3);
t1136 = legFrame(2,3) + qJ(1,2);
t1137 = legFrame(1,3) + qJ(1,1);
t1216 = 0.2e1 * pkin(1);
t1153 = xDP(2);
t1154 = xDP(1);
t1124 = pkin(7) + t1137;
t1121 = qJ(3,1) + t1124;
t1112 = cos(t1121);
t1169 = -pkin(1) * cos(t1137) - pkin(2) * cos(t1124) - pkin(3) * t1112;
t1109 = sin(t1121);
t1172 = -pkin(1) * sin(t1137) - pkin(2) * sin(t1124) - pkin(3) * t1109;
t1061 = t1172 * t1153 + t1169 * t1154;
t1156 = 0.1e1 / pkin(3);
t1058 = t1061 * t1156 * t1098;
t1076 = t1109 * t1153 + t1112 * t1154;
t1064 = t1076 * t1098;
t1049 = t1064 + t1058 / 0.2e1;
t1055 = t1064 + t1058;
t1134 = cos(t1141);
t1138 = pkin(2) ^ 2 + t1157;
t1155 = pkin(3) ^ 2;
t1215 = (t1049 * t1152 * t1218 + t1138 * t1064 + t1055 * t1155 + (pkin(3) * t1049 * t1134 + t1064 * t1222) * t1216) * t1156;
t1122 = pkin(7) + t1135;
t1119 = qJ(3,3) + t1122;
t1110 = cos(t1119);
t1171 = -pkin(1) * cos(t1135) - pkin(2) * cos(t1122) - pkin(3) * t1110;
t1107 = sin(t1119);
t1174 = -pkin(1) * sin(t1135) - pkin(2) * sin(t1122) - pkin(3) * t1107;
t1060 = t1174 * t1153 + t1171 * t1154;
t1057 = t1060 * t1156 * t1092;
t1074 = t1107 * t1153 + t1110 * t1154;
t1062 = t1074 * t1092;
t1047 = t1062 + t1057 / 0.2e1;
t1051 = t1057 + t1062;
t1132 = cos(t1139);
t1214 = (t1047 * t1150 * t1218 + t1138 * t1062 + t1051 * t1155 + (pkin(3) * t1047 * t1132 + t1062 * t1222) * t1216) * t1156;
t1123 = pkin(7) + t1136;
t1120 = qJ(3,2) + t1123;
t1111 = cos(t1120);
t1170 = -pkin(1) * cos(t1136) - pkin(2) * cos(t1123) - pkin(3) * t1111;
t1108 = sin(t1120);
t1173 = -pkin(1) * sin(t1136) - pkin(2) * sin(t1123) - pkin(3) * t1108;
t1059 = t1173 * t1153 + t1170 * t1154;
t1056 = t1059 * t1156 * t1095;
t1075 = t1108 * t1153 + t1111 * t1154;
t1063 = t1075 * t1095;
t1048 = t1063 + t1056 / 0.2e1;
t1050 = t1056 + t1063;
t1133 = cos(t1140);
t1213 = (t1048 * t1151 * t1218 + t1138 * t1063 + t1050 * t1155 + (pkin(3) * t1048 * t1133 + t1063 * t1222) * t1216) * t1156;
t1212 = t1059 * t1096;
t1211 = t1060 * t1093;
t1210 = t1061 * t1099;
t1209 = t1074 ^ 2 * t1092 * t1093;
t1208 = t1075 ^ 2 * t1095 * t1096;
t1207 = t1076 ^ 2 * t1098 * t1099;
t1206 = t1074 * t1093;
t1205 = t1075 * t1096;
t1204 = t1076 * t1099;
t1131 = pkin(1) * t1143 + pkin(2);
t1089 = t1131 * t1147 + t1150 * t1223;
t1203 = t1089 * t1092;
t1090 = t1131 * t1148 + t1151 * t1223;
t1202 = t1090 * t1095;
t1091 = t1131 * t1149 + t1152 * t1223;
t1201 = t1091 * t1098;
t1200 = t1092 * t1107;
t1199 = t1092 * t1110;
t1198 = t1095 * t1108;
t1197 = t1095 * t1111;
t1196 = t1098 * t1109;
t1195 = t1098 * t1112;
t1194 = t1142 * t1147;
t1193 = t1142 * t1148;
t1192 = t1142 * t1149;
t1191 = 2 * MDP(6);
t1190 = 2 * MDP(7);
t1041 = -pkin(3) * t1051 + (-pkin(1) * t1132 - t1221) * t1062;
t1045 = t1051 * t1211;
t1086 = -pkin(1) * t1194 + t1150 * t1131;
t1176 = t1051 * (pkin(3) + t1086) / t1089 * t1057;
t1189 = (t1045 - t1176 / 0.2e1 + (-t1041 - t1214 / 0.2e1) * t1206) * t1203;
t1077 = -t1221 + (-t1143 * t1150 + t1194) * pkin(1);
t1188 = t1092 * (-t1176 + 0.2e1 * t1045 + (-0.2e1 * t1041 - t1214) * t1206) * t1077;
t1043 = -pkin(3) * t1050 + (-pkin(1) * t1133 - t1220) * t1063;
t1044 = t1050 * t1212;
t1087 = -pkin(1) * t1193 + t1151 * t1131;
t1177 = t1050 * (pkin(3) + t1087) / t1090 * t1056;
t1187 = (t1044 - t1177 / 0.2e1 + (-t1043 - t1213 / 0.2e1) * t1205) * t1202;
t1078 = -t1220 + (-t1143 * t1151 + t1193) * pkin(1);
t1186 = t1095 * (-t1177 + 0.2e1 * t1044 + (-0.2e1 * t1043 - t1213) * t1205) * t1078;
t1042 = -pkin(3) * t1055 + (-pkin(1) * t1134 - t1219) * t1064;
t1046 = t1055 * t1210;
t1088 = -pkin(1) * t1192 + t1152 * t1131;
t1175 = t1055 * (pkin(3) + t1088) / t1091 * t1058;
t1185 = (t1046 - t1175 / 0.2e1 + (-t1042 - t1215 / 0.2e1) * t1204) * t1201;
t1079 = -t1219 + (-t1143 * t1152 + t1192) * pkin(1);
t1184 = t1098 * (-t1175 + 0.2e1 * t1046 + (-0.2e1 * t1042 - t1215) * t1204) * t1079;
t1183 = t1047 * t1089 * t1211;
t1182 = (0.2e1 * t1062 + t1057) * t1077 * t1211;
t1181 = t1048 * t1090 * t1212;
t1180 = (0.2e1 * t1063 + t1056) * t1078 * t1212;
t1179 = t1049 * t1091 * t1210;
t1178 = (0.2e1 * t1064 + t1058) * t1079 * t1210;
t1032 = -t1176 + t1045 + (-t1041 - t1214) * t1206;
t1038 = -t1041 * t1206 + t1045;
t1166 = t1032 * t1092 * MDP(5) + (t1092 * t1086 * t1038 + t1089 * t1209) * MDP(6) + (-t1038 * t1203 + t1086 * t1209) * MDP(7);
t1033 = -t1175 + t1046 + (-t1042 - t1215) * t1204;
t1039 = -t1042 * t1204 + t1046;
t1165 = t1033 * t1098 * MDP(5) + (t1098 * t1088 * t1039 + t1091 * t1207) * MDP(6) + (-t1039 * t1201 + t1088 * t1207) * MDP(7);
t1034 = -t1177 + t1044 + (-t1043 - t1213) * t1205;
t1040 = -t1043 * t1205 + t1044;
t1164 = t1034 * t1095 * MDP(5) + (t1095 * t1087 * t1040 + t1090 * t1208) * MDP(6) + (-t1040 * t1202 + t1087 * t1208) * MDP(7);
t1 = [(t1032 * t1199 + t1033 * t1195 + t1034 * t1197) * MDP(5) + (-t1110 * t1188 - t1111 * t1186 - t1112 * t1184) * MDP(6) + (-t1110 * t1189 - t1111 * t1187 - t1112 * t1185) * t1190 + ((-t1110 * t1183 - t1111 * t1181 - t1112 * t1179) * t1191 + (t1110 * t1182 + t1111 * t1180 + t1112 * t1178) * MDP(7) + t1165 * t1169 + t1164 * t1170 + t1166 * t1171) * t1156 + t1224 * (t1038 * t1199 + t1039 * t1195 + t1040 * t1197); (t1032 * t1200 + t1033 * t1196 + t1034 * t1198) * MDP(5) + (-t1107 * t1188 - t1108 * t1186 - t1109 * t1184) * MDP(6) + (-t1107 * t1189 - t1108 * t1187 - t1109 * t1185) * t1190 + ((-t1107 * t1183 - t1108 * t1181 - t1109 * t1179) * t1191 + (t1107 * t1182 + t1108 * t1180 + t1109 * t1178) * MDP(7) + t1165 * t1172 + t1164 * t1173 + t1166 * t1174) * t1156 + t1224 * (t1038 * t1200 + t1039 * t1196 + t1040 * t1198); 0;];
taucX  = t1;
