% Calculate Gravitation load for parallel robot
% P3RRPRR8V2G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 21:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V2G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:12:10
% EndTime: 2020-08-06 21:12:11
% DurationCPUTime: 0.70s
% Computational Cost: add. (771->168), mult. (1143->276), div. (36->9), fcn. (756->23), ass. (0->127)
t1099 = cos(qJ(2,3));
t1174 = 0.2e1 * t1099 ^ 2;
t1101 = cos(qJ(2,2));
t1173 = 0.2e1 * t1101 ^ 2;
t1103 = cos(qJ(2,1));
t1172 = 0.2e1 * t1103 ^ 2;
t1089 = -qJ(3,1) - pkin(5);
t1088 = -qJ(3,2) - pkin(5);
t1087 = -qJ(3,3) - pkin(5);
t1094 = sin(qJ(1,3));
t1078 = pkin(6) - t1087;
t1100 = cos(qJ(1,3));
t1121 = pkin(1) * t1094 - t1100 * t1078;
t1086 = cos(pkin(7));
t1081 = t1086 ^ 2;
t1141 = pkin(3) * (t1081 - 0.1e1);
t1085 = sin(pkin(7));
t1093 = sin(qJ(2,3));
t1145 = t1085 * t1093;
t1171 = pkin(3) * (t1094 * t1141 + t1121 * t1145);
t1096 = sin(qJ(1,2));
t1079 = pkin(6) - t1088;
t1102 = cos(qJ(1,2));
t1120 = pkin(1) * t1096 - t1102 * t1079;
t1095 = sin(qJ(2,2));
t1144 = t1085 * t1095;
t1170 = pkin(3) * (t1096 * t1141 + t1120 * t1144);
t1098 = sin(qJ(1,1));
t1080 = pkin(6) - t1089;
t1104 = cos(qJ(1,1));
t1119 = pkin(1) * t1098 - t1104 * t1080;
t1097 = sin(qJ(2,1));
t1143 = t1085 * t1097;
t1169 = pkin(3) * (t1098 * t1141 + t1119 * t1143);
t1168 = t1085 * pkin(3);
t1167 = t1086 * pkin(3);
t1090 = legFrame(3,2);
t1065 = sin(t1090);
t1068 = cos(t1090);
t1044 = t1068 * g(1) - t1065 * g(2);
t1062 = 0.1e1 / t1078;
t1166 = (t1094 * g(3) - t1100 * t1044) * t1062;
t1091 = legFrame(2,2);
t1066 = sin(t1091);
t1069 = cos(t1091);
t1045 = t1069 * g(1) - t1066 * g(2);
t1063 = 0.1e1 / t1079;
t1165 = (t1096 * g(3) - t1102 * t1045) * t1063;
t1092 = legFrame(1,2);
t1067 = sin(t1092);
t1070 = cos(t1092);
t1046 = t1070 * g(1) - t1067 * g(2);
t1064 = 0.1e1 / t1080;
t1164 = (t1098 * g(3) - t1104 * t1046) * t1064;
t1040 = -m(3) * pkin(2) - mrSges(3,1) * t1086 + mrSges(3,2) * t1085 - mrSges(2,1);
t1041 = t1065 * g(1) + t1068 * g(2);
t1050 = t1085 * mrSges(3,1) + t1086 * mrSges(3,2) + mrSges(2,2);
t1118 = -g(3) * t1100 - t1044 * t1094;
t1125 = pkin(3) * cos(qJ(2,3) + pkin(7)) + t1099 * pkin(2);
t1163 = 0.1e1 / t1125 * ((t1118 * t1040 + t1041 * t1050) * t1093 - t1099 * (-t1041 * t1040 + t1118 * t1050));
t1042 = t1066 * g(1) + t1069 * g(2);
t1117 = -g(3) * t1102 - t1045 * t1096;
t1124 = pkin(3) * cos(qJ(2,2) + pkin(7)) + t1101 * pkin(2);
t1162 = 0.1e1 / t1124 * ((t1117 * t1040 + t1042 * t1050) * t1095 - t1101 * (-t1042 * t1040 + t1117 * t1050));
t1043 = t1067 * g(1) + t1070 * g(2);
t1116 = -g(3) * t1104 - t1046 * t1098;
t1123 = pkin(3) * cos(qJ(2,1) + pkin(7)) + t1103 * pkin(2);
t1161 = 0.1e1 / t1123 * ((t1116 * t1040 + t1043 * t1050) * t1097 - t1103 * (-t1043 * t1040 + t1116 * t1050));
t1056 = pkin(2) + t1167;
t1160 = t1056 * t1068;
t1159 = t1056 * t1069;
t1158 = t1056 * t1070;
t1055 = mrSges(1,1) + (m(2) + m(3)) * pkin(1);
t1054 = t1055 * g(3);
t1122 = -m(2) * pkin(5) + mrSges(1,2) - mrSges(2,3) - mrSges(3,3);
t1112 = t1087 * m(3) + t1122;
t1115 = -t1040 * t1099 - t1050 * t1093;
t1157 = t1062 * (t1054 * t1094 + (t1115 * t1094 + t1112 * t1100) * g(3) + ((-t1055 - t1115) * t1100 + t1112 * t1094) * t1044);
t1111 = t1088 * m(3) + t1122;
t1114 = -t1040 * t1101 - t1050 * t1095;
t1156 = t1063 * (t1054 * t1096 + (t1114 * t1096 + t1111 * t1102) * g(3) + ((-t1055 - t1114) * t1102 + t1111 * t1096) * t1045);
t1110 = t1089 * m(3) + t1122;
t1113 = -t1040 * t1103 - t1050 * t1097;
t1155 = t1064 * (t1054 * t1098 + (t1113 * t1098 + t1110 * t1104) * g(3) + ((-t1055 - t1113) * t1104 + t1110 * t1098) * t1046);
t1154 = t1065 * t1056;
t1153 = t1065 * t1094;
t1152 = t1066 * t1056;
t1151 = t1066 * t1096;
t1150 = t1067 * t1056;
t1149 = t1067 * t1098;
t1148 = t1068 * t1094;
t1147 = t1069 * t1096;
t1146 = t1070 * t1098;
t1142 = pkin(2) * t1167;
t1140 = t1068 * t1168;
t1139 = t1069 * t1168;
t1138 = t1070 * t1168;
t1137 = t1065 * t1168;
t1136 = t1066 * t1168;
t1135 = t1067 * t1168;
t1134 = pkin(3) * t1145;
t1133 = pkin(3) * t1144;
t1132 = pkin(3) * t1143;
t1036 = 0.1e1 / (t1056 * t1099 - t1134);
t1131 = t1036 * t1157;
t1037 = 0.1e1 / (t1056 * t1101 - t1133);
t1130 = t1037 * t1156;
t1038 = 0.1e1 / (t1056 * t1103 - t1132);
t1129 = t1038 * t1155;
t1128 = t1036 * t1166;
t1127 = t1037 * t1165;
t1126 = t1038 * t1164;
t1107 = pkin(3) ^ 2;
t1108 = pkin(2) ^ 2;
t1109 = 0.2e1 * t1081 * t1107 - t1107 + t1108 + 0.2e1 * t1142;
t1057 = pkin(1) * t1168;
t1053 = pkin(1) * t1097 - t1168;
t1052 = pkin(1) * t1095 - t1168;
t1051 = pkin(1) * t1093 - t1168;
t1039 = t1142 + t1108 / 0.2e1 + (t1081 - 0.1e1 / 0.2e1) * t1107;
t1032 = -0.2e1 * t1098 * t1132 + t1119;
t1031 = -0.2e1 * t1096 * t1133 + t1120;
t1030 = -0.2e1 * t1094 * t1134 + t1121;
t1029 = t1109 * t1097 + t1057;
t1028 = t1109 * t1095 + t1057;
t1027 = t1109 * t1093 + t1057;
t1 = [((t1056 * t1146 + t1135) * t1103 + t1097 * (-t1098 * t1138 + t1150)) * t1129 + t1067 * t1161 + ((t1056 * t1147 + t1136) * t1101 + t1095 * (-t1096 * t1139 + t1152)) * t1130 + t1066 * t1162 + ((t1056 * t1148 + t1137) * t1099 + t1093 * (-t1094 * t1140 + t1154)) * t1131 + t1065 * t1163 - g(1) * m(4) + (-((t1039 * t1146 + t1056 * t1135) * t1172 + (t1067 * t1029 + t1032 * t1158) * t1103 - t1070 * t1169 + t1053 * t1150) * t1126 - ((t1039 * t1147 + t1056 * t1136) * t1173 + (t1066 * t1028 + t1031 * t1159) * t1101 - t1069 * t1170 + t1052 * t1152) * t1127 - ((t1039 * t1148 + t1056 * t1137) * t1174 + (t1065 * t1027 + t1030 * t1160) * t1099 - t1068 * t1171 + t1051 * t1154) * t1128) * m(3); ((-t1056 * t1149 + t1138) * t1103 + (t1098 * t1135 + t1158) * t1097) * t1129 + t1070 * t1161 + ((-t1056 * t1151 + t1139) * t1101 + (t1096 * t1136 + t1159) * t1095) * t1130 + t1069 * t1162 + ((-t1056 * t1153 + t1140) * t1099 + (t1094 * t1137 + t1160) * t1093) * t1131 + t1068 * t1163 - g(2) * m(4) + (-((-t1039 * t1149 + t1056 * t1138) * t1172 + (t1029 * t1070 - t1032 * t1150) * t1103 + t1067 * t1169 + t1053 * t1158) * t1126 - ((-t1039 * t1151 + t1056 * t1139) * t1173 + (t1028 * t1069 - t1031 * t1152) * t1101 + t1066 * t1170 + t1052 * t1159) * t1127 - ((-t1039 * t1153 + t1056 * t1140) * t1174 + (t1027 * t1068 - t1030 * t1154) * t1099 + t1065 * t1171 + t1051 * t1160) * t1128) * m(3); -g(3) * m(4) + t1100 * t1157 + t1102 * t1156 + t1104 * t1155 + (-(t1098 * t1080 + (pkin(1) + t1123) * t1104) * t1164 - (t1096 * t1079 + (pkin(1) + t1124) * t1102) * t1165 - (t1094 * t1078 + (pkin(1) + t1125) * t1100) * t1166) * m(3);];
taugX  = t1;
