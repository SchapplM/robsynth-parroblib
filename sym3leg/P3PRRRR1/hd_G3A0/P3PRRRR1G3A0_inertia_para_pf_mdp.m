% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRRR1G3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR1G3A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3PRRRR1G3A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3A0_inertia_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR1G3A0_inertia_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:31
% EndTime: 2020-03-09 21:02:33
% DurationCPUTime: 1.98s
% Computational Cost: add. (519->232), mult. (1769->430), div. (1002->17), fcn. (2034->18), ass. (0->195)
t1062 = legFrame(3,2);
t1032 = cos(t1062);
t1063 = legFrame(2,2);
t1033 = cos(t1063);
t1064 = legFrame(1,2);
t1034 = cos(t1064);
t1076 = cos(qJ(2,1));
t1069 = sin(qJ(3,1));
t1044 = t1069 ^ 2;
t1070 = sin(qJ(2,1));
t1046 = 0.1e1 / t1070 ^ 2;
t1075 = cos(qJ(3,1));
t1091 = t1075 ^ 2;
t1058 = 0.1e1 / t1091;
t1185 = t1046 * t1058;
t1150 = t1044 * t1185;
t1116 = t1076 * t1150;
t1074 = cos(qJ(2,2));
t1067 = sin(qJ(3,2));
t1040 = t1067 ^ 2;
t1068 = sin(qJ(2,2));
t1042 = 0.1e1 / t1068 ^ 2;
t1073 = cos(qJ(3,2));
t1088 = t1073 ^ 2;
t1053 = 0.1e1 / t1088;
t1190 = t1042 * t1053;
t1158 = t1040 * t1190;
t1121 = t1074 * t1158;
t1072 = cos(qJ(2,3));
t1065 = sin(qJ(3,3));
t1036 = t1065 ^ 2;
t1066 = sin(qJ(2,3));
t1038 = 0.1e1 / t1066 ^ 2;
t1071 = cos(qJ(3,3));
t1085 = t1071 ^ 2;
t1048 = 0.1e1 / t1085;
t1195 = t1038 * t1048;
t1166 = t1036 * t1195;
t1126 = t1072 * t1166;
t1221 = t1032 * t1126 + t1033 * t1121 + t1034 * t1116;
t1029 = sin(t1062);
t1220 = t1029 * t1032;
t1030 = sin(t1063);
t1219 = t1030 * t1033;
t1031 = sin(t1064);
t1218 = t1031 * t1034;
t1217 = 2 * MDP(6);
t1047 = 0.1e1 / t1071;
t1052 = 0.1e1 / t1073;
t1057 = 0.1e1 / t1075;
t1077 = 1 / pkin(2);
t1216 = 2 * t1077;
t1078 = 1 / (pkin(2) ^ 2);
t1215 = MDP(2) * t1078;
t1214 = MDP(3) * t1077;
t1213 = MDP(4) * t1077;
t1212 = MDP(5) * t1078;
t1211 = MDP(7) * t1078;
t1210 = MDP(8) * t1078;
t1203 = t1029 * t1072;
t1017 = t1032 * t1066 - t1203;
t1037 = 0.1e1 / t1066;
t1209 = t1017 * t1037;
t1200 = t1032 * t1072;
t1018 = t1066 * t1029 + t1200;
t1208 = t1018 * t1037;
t1202 = t1030 * t1074;
t1019 = t1033 * t1068 - t1202;
t1041 = 0.1e1 / t1068;
t1207 = t1019 * t1041;
t1199 = t1033 * t1074;
t1020 = t1068 * t1030 + t1199;
t1206 = t1020 * t1041;
t1201 = t1031 * t1076;
t1021 = t1034 * t1070 - t1201;
t1045 = 0.1e1 / t1070;
t1205 = t1021 * t1045;
t1198 = t1034 * t1076;
t1022 = t1070 * t1031 + t1198;
t1204 = t1022 * t1045;
t1197 = t1036 * t1048;
t1196 = t1037 * t1047;
t1051 = t1072 ^ 2;
t1194 = t1038 * t1051;
t1193 = t1038 * t1072;
t1192 = t1040 * t1053;
t1191 = t1041 * t1052;
t1056 = t1074 ^ 2;
t1189 = t1042 * t1056;
t1188 = t1042 * t1074;
t1187 = t1044 * t1058;
t1186 = t1045 * t1057;
t1061 = t1076 ^ 2;
t1184 = t1046 * t1061;
t1183 = t1046 * t1076;
t1182 = t1047 * t1065;
t1181 = t1048 * t1065;
t1049 = t1047 * t1048;
t1180 = t1049 * t1072;
t1179 = t1052 * t1067;
t1178 = t1053 * t1067;
t1054 = t1052 * t1053;
t1177 = t1054 * t1074;
t1176 = t1057 * t1069;
t1175 = t1058 * t1069;
t1059 = t1057 * t1058;
t1174 = t1059 * t1076;
t1173 = t1078 * t1217;
t1172 = t1029 * t1196;
t1171 = t1030 * t1191;
t1170 = t1031 * t1186;
t1169 = t1032 * t1196;
t1168 = t1033 * t1191;
t1167 = t1034 * t1186;
t1165 = t1051 * t1197;
t1164 = t1037 * t1182;
t1163 = t1037 * t1181;
t1162 = t1038 * t1182;
t1161 = t1049 * t1194;
t1160 = 0.1e1 / t1085 ^ 2 * t1194;
t1159 = t1065 * t1193;
t1157 = t1056 * t1192;
t1156 = t1041 * t1179;
t1155 = t1041 * t1178;
t1154 = t1042 * t1179;
t1153 = t1054 * t1189;
t1152 = 0.1e1 / t1088 ^ 2 * t1189;
t1151 = t1067 * t1188;
t1149 = t1061 * t1187;
t1148 = t1045 * t1176;
t1147 = t1045 * t1175;
t1146 = t1046 * t1176;
t1145 = t1059 * t1184;
t1144 = 0.1e1 / t1091 ^ 2 * t1184;
t1143 = t1069 * t1183;
t1115 = t1057 * t1143;
t1111 = t1031 * t1115;
t1120 = t1052 * t1151;
t1112 = t1030 * t1120;
t1125 = t1047 * t1159;
t1113 = t1029 * t1125;
t1142 = (t1111 + t1112 + t1113) * t1077;
t1141 = t1221 * t1077;
t1140 = t1017 * t1029 * t1193;
t1139 = t1017 * t1162;
t1138 = t1018 * t1032 * t1193;
t1137 = t1019 * t1030 * t1188;
t1136 = t1019 * t1154;
t1135 = t1020 * t1033 * t1188;
t1134 = t1021 * t1031 * t1183;
t1133 = t1021 * t1146;
t1132 = t1022 * t1034 * t1183;
t1131 = t1195 * t1220;
t1130 = t1190 * t1219;
t1129 = t1185 * t1218;
t1035 = t1065 * t1036;
t1128 = t1035 * t1038 * t1180;
t1127 = t1038 * t1165;
t1124 = t1049 * t1159;
t1039 = t1067 * t1040;
t1123 = t1039 * t1042 * t1177;
t1122 = t1042 * t1157;
t1119 = t1054 * t1151;
t1043 = t1069 * t1044;
t1118 = t1043 * t1046 * t1174;
t1117 = t1046 * t1149;
t1114 = t1059 * t1143;
t1109 = t1032 * t1125;
t1107 = t1033 * t1120;
t1105 = t1034 * t1115;
t1104 = t1017 * t1032 - t1018 * t1029;
t1103 = t1019 * t1033 - t1020 * t1030;
t1102 = t1021 * t1034 - t1022 * t1031;
t1101 = (-t1018 * t1051 - t1200) * t1038;
t1100 = (-t1020 * t1056 - t1199) * t1042;
t1099 = (-t1022 * t1061 - t1198) * t1046;
t1094 = -t1029 * t1126 - t1030 * t1121 - t1031 * t1116;
t1098 = (t1029 * t1163 + t1030 * t1155 + t1031 * t1147) * t1211 + (t1170 + t1171 + t1172) * t1210 + ((-t1021 * t1061 + t1201) * t1046 * t1175 + (-t1019 * t1056 + t1202) * t1042 * t1178 + (-t1017 * t1051 + t1203) * t1038 * t1181) * t1214 + ((t1021 * t1076 - t1031) * t1147 + (t1019 * t1074 - t1030) * t1155 + (t1017 * t1072 - t1029) * t1163) * t1213 + (t1133 + t1136 + t1139) * MDP(1) + t1094 * t1173 + (-t1029 * t1128 - t1030 * t1123 - t1031 * t1118) * t1212 + (-t1029 * t1124 - t1030 * t1119 - t1031 * t1114) * t1215;
t1097 = (t1032 * t1124 + t1033 * t1119 + t1034 * t1114) * t1215 + (-t1032 * t1163 - t1033 * t1155 - t1034 * t1147) * t1211 + (-t1167 - t1168 - t1169) * t1210 + (t1099 * t1175 + t1100 * t1178 + t1101 * t1181) * t1214 + ((t1022 * t1076 + t1034) * t1147 + (t1020 * t1074 + t1033) * t1155 + (t1018 * t1072 + t1032) * t1163) * t1213 + (t1018 * t1162 + t1020 * t1154 + t1022 * t1146) * MDP(1) + t1221 * t1173 + (t1032 * t1128 + t1033 * t1123 + t1034 * t1118) * t1212;
t1096 = t1035 * t1161 + t1039 * t1153 + t1043 * t1145;
t1095 = t1036 * t1037 * t1180 + t1040 * t1041 * t1177 + t1044 * t1045 * t1174;
t1028 = t1034 ^ 2;
t1027 = t1033 ^ 2;
t1026 = t1032 ^ 2;
t1025 = t1031 ^ 2;
t1024 = t1030 ^ 2;
t1023 = t1029 ^ 2;
t1010 = (t1045 * t1149 - t1070) * t1077;
t1009 = (t1041 * t1157 - t1068) * t1077;
t1008 = (t1037 * t1165 - t1066) * t1077;
t1007 = (-t1045 * t1061 - t1070) * t1077 * t1176;
t1006 = (-t1041 * t1056 - t1068) * t1077 * t1179;
t1005 = (-t1037 * t1051 - t1066) * t1077 * t1182;
t988 = (t1018 * t1038 * t1017 + t1020 * t1042 * t1019 + t1022 * t1046 * t1021) * MDP(1) + ((-t1129 - t1130 - t1131) * MDP(2) + (-t1036 * t1131 - t1040 * t1130 - t1044 * t1129) * MDP(5) + (-t1146 * t1218 - t1154 * t1219 - t1162 * t1220) * t1217) * t1078 + ((t1102 * t1186 + t1103 * t1191 + t1104 * t1196) * MDP(4) + (-MDP(10) + (MDP(11) * t1069 - MDP(3)) * t1057) * t1102 * t1183 + (-MDP(10) + (MDP(11) * t1067 - MDP(3)) * t1052) * t1103 * t1188 + (-MDP(10) + (MDP(11) * t1065 - MDP(3)) * t1047) * t1104 * t1193) * t1077;
t1 = [(t1018 ^ 2 * t1038 + t1020 ^ 2 * t1042 + t1022 ^ 2 * t1046) * MDP(1) + MDP(12) + ((t1026 * t1195 + t1027 * t1190 + t1028 * t1185) * MDP(2) + (t1026 * t1166 + t1027 * t1158 + t1028 * t1150) * MDP(5) + (t1026 * t1162 + t1027 * t1154 + t1028 * t1146) * t1217) * t1078 + ((-t1047 * t1138 - t1052 * t1135 - t1057 * t1132) * MDP(3) + (t1018 * t1169 + t1020 * t1168 + t1022 * t1167) * MDP(4) + (-t1132 - t1135 - t1138) * MDP(10) + (t1018 * t1109 + t1020 * t1107 + t1022 * t1105) * MDP(11)) * t1216; t988; (t1008 * t1208 + t1009 * t1206 + t1010 * t1204 + t1141) * MDP(11) + t1097 + (t1005 * t1208 + t1006 * t1206 + t1007 * t1204 + (-t1105 - t1107 - t1109) * t1077) * MDP(10); t988; (t1017 ^ 2 * t1038 + t1019 ^ 2 * t1042 + t1021 ^ 2 * t1046) * MDP(1) + MDP(12) + ((t1023 * t1195 + t1024 * t1190 + t1025 * t1185) * MDP(2) + (t1023 * t1166 + t1024 * t1158 + t1025 * t1150) * MDP(5) + (t1023 * t1162 + t1024 * t1154 + t1025 * t1146) * t1217) * t1078 + ((t1047 * t1140 + t1052 * t1137 + t1057 * t1134) * MDP(3) + (-t1017 * t1172 - t1019 * t1171 - t1021 * t1170) * MDP(4) + (t1134 + t1137 + t1140) * MDP(10) + (-t1017 * t1113 - t1019 * t1112 - t1021 * t1111) * MDP(11)) * t1216; (t1005 * t1209 + t1006 * t1207 + t1007 * t1205 + t1142) * MDP(10) + t1098 + (t1008 * t1209 + t1009 * t1207 + t1010 * t1205 + t1094 * t1077) * MDP(11); t1097 + (t1141 + (t1018 * t1127 + t1020 * t1122 + t1022 * t1117 - t1018 - t1020 - t1022) * t1077) * MDP(11) + ((-t1022 + t1099) * t1176 + (-t1020 + t1100) * t1179 + (-t1018 + t1101) * t1182) * t1077 * MDP(10); t1142 * MDP(10) + ((-t1017 * t1182 - t1019 * t1179 - t1021 * t1176 - t1051 * t1139 - t1056 * t1136 - t1061 * t1133) * MDP(10) + (t1017 * t1127 + t1019 * t1122 + t1021 * t1117 - t1017 - t1019 - t1021 + t1094) * MDP(11)) * t1077 + t1098; (t1150 + t1158 + t1166) * MDP(1) + (t1005 * t1164 + t1006 * t1156 + t1007 * t1148) * MDP(10) + (t1008 * t1164 + t1009 * t1156 + t1010 * t1148) * MDP(11) + MDP(12) + ((-t1117 - t1122 - t1127 - t1187 - t1192 - t1197) * MDP(10) + (t1096 - t1176 - t1179 - t1182) * MDP(11) + 0.2e1 * (-t1036 * t1161 - t1040 * t1153 - t1044 * t1145) * MDP(3) + 0.2e1 * t1095 * MDP(4)) * t1077 + ((t1036 * t1160 + t1040 * t1152 + t1044 * t1144) * MDP(2) + (t1036 ^ 2 * t1160 + t1040 ^ 2 * t1152 + t1044 ^ 2 * t1144) * MDP(5) + (t1048 + t1053 + t1058) * MDP(9) + t1096 * t1217 - 0.2e1 * t1095 * MDP(7) + 0.2e1 * (-t1072 * t1163 - t1074 * t1155 - t1076 * t1147) * MDP(8)) * t1078;];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
