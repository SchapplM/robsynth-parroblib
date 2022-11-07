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
% Datum: 2022-11-07 13:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
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
% StartTime: 2022-11-07 13:06:07
% EndTime: 2022-11-07 13:06:08
% DurationCPUTime: 0.65s
% Computational Cost: add. (771->168), mult. (1143->276), div. (36->9), fcn. (756->23), ass. (0->127)
t1091 = cos(qJ(2,3));
t1166 = 0.2e1 * t1091 ^ 2;
t1093 = cos(qJ(2,2));
t1165 = 0.2e1 * t1093 ^ 2;
t1095 = cos(qJ(2,1));
t1164 = 0.2e1 * t1095 ^ 2;
t1081 = qJ(3,1) + pkin(5);
t1080 = qJ(3,2) + pkin(5);
t1079 = qJ(3,3) + pkin(5);
t1086 = sin(qJ(1,3));
t1070 = pkin(6) + t1079;
t1092 = cos(qJ(1,3));
t1113 = pkin(1) * t1086 - t1092 * t1070;
t1078 = cos(pkin(7));
t1073 = t1078 ^ 2;
t1133 = pkin(3) * (t1073 - 0.1e1);
t1077 = sin(pkin(7));
t1085 = sin(qJ(2,3));
t1137 = t1077 * t1085;
t1163 = (t1086 * t1133 + t1113 * t1137) * pkin(3);
t1088 = sin(qJ(1,2));
t1071 = pkin(6) + t1080;
t1094 = cos(qJ(1,2));
t1112 = pkin(1) * t1088 - t1094 * t1071;
t1087 = sin(qJ(2,2));
t1136 = t1077 * t1087;
t1162 = (t1088 * t1133 + t1112 * t1136) * pkin(3);
t1090 = sin(qJ(1,1));
t1072 = pkin(6) + t1081;
t1096 = cos(qJ(1,1));
t1111 = pkin(1) * t1090 - t1096 * t1072;
t1089 = sin(qJ(2,1));
t1135 = t1077 * t1089;
t1161 = (t1090 * t1133 + t1111 * t1135) * pkin(3);
t1160 = t1077 * pkin(3);
t1159 = t1078 * pkin(3);
t1082 = legFrame(3,2);
t1057 = sin(t1082);
t1060 = cos(t1082);
t1036 = t1060 * g(1) - t1057 * g(2);
t1054 = 0.1e1 / t1070;
t1158 = (t1086 * g(3) - t1092 * t1036) * t1054;
t1083 = legFrame(2,2);
t1058 = sin(t1083);
t1061 = cos(t1083);
t1037 = t1061 * g(1) - t1058 * g(2);
t1055 = 0.1e1 / t1071;
t1157 = (t1088 * g(3) - t1094 * t1037) * t1055;
t1084 = legFrame(1,2);
t1059 = sin(t1084);
t1062 = cos(t1084);
t1038 = t1062 * g(1) - t1059 * g(2);
t1056 = 0.1e1 / t1072;
t1156 = (t1090 * g(3) - t1096 * t1038) * t1056;
t1032 = -m(3) * pkin(2) - mrSges(3,1) * t1078 + mrSges(3,2) * t1077 - mrSges(2,1);
t1033 = t1057 * g(1) + t1060 * g(2);
t1042 = t1077 * mrSges(3,1) + t1078 * mrSges(3,2) + mrSges(2,2);
t1110 = -g(3) * t1092 - t1036 * t1086;
t1117 = pkin(3) * cos(qJ(2,3) + pkin(7)) + t1091 * pkin(2);
t1155 = 0.1e1 / t1117 * ((t1110 * t1032 + t1033 * t1042) * t1085 - (-t1033 * t1032 + t1110 * t1042) * t1091);
t1034 = t1058 * g(1) + t1061 * g(2);
t1109 = -g(3) * t1094 - t1037 * t1088;
t1116 = pkin(3) * cos(qJ(2,2) + pkin(7)) + t1093 * pkin(2);
t1154 = 0.1e1 / t1116 * ((t1109 * t1032 + t1034 * t1042) * t1087 - (-t1034 * t1032 + t1109 * t1042) * t1093);
t1035 = t1059 * g(1) + t1062 * g(2);
t1108 = -g(3) * t1096 - t1038 * t1090;
t1115 = pkin(3) * cos(qJ(2,1) + pkin(7)) + t1095 * pkin(2);
t1153 = 0.1e1 / t1115 * ((t1108 * t1032 + t1035 * t1042) * t1089 - (-t1035 * t1032 + t1108 * t1042) * t1095);
t1048 = pkin(2) + t1159;
t1152 = t1048 * t1060;
t1151 = t1048 * t1061;
t1150 = t1048 * t1062;
t1047 = mrSges(1,1) + (m(2) + m(3)) * pkin(1);
t1046 = t1047 * g(3);
t1114 = -m(2) * pkin(5) + mrSges(1,2) - mrSges(2,3) - mrSges(3,3);
t1104 = -t1079 * m(3) + t1114;
t1107 = -t1032 * t1091 - t1042 * t1085;
t1149 = t1054 * (t1046 * t1086 + (t1107 * t1086 + t1104 * t1092) * g(3) + ((-t1047 - t1107) * t1092 + t1104 * t1086) * t1036);
t1103 = -t1080 * m(3) + t1114;
t1106 = -t1032 * t1093 - t1042 * t1087;
t1148 = t1055 * (t1046 * t1088 + (t1106 * t1088 + t1103 * t1094) * g(3) + ((-t1047 - t1106) * t1094 + t1103 * t1088) * t1037);
t1102 = -t1081 * m(3) + t1114;
t1105 = -t1032 * t1095 - t1042 * t1089;
t1147 = t1056 * (t1046 * t1090 + (t1105 * t1090 + t1102 * t1096) * g(3) + ((-t1047 - t1105) * t1096 + t1102 * t1090) * t1038);
t1146 = t1057 * t1048;
t1145 = t1057 * t1086;
t1144 = t1058 * t1048;
t1143 = t1058 * t1088;
t1142 = t1059 * t1048;
t1141 = t1059 * t1090;
t1140 = t1060 * t1086;
t1139 = t1061 * t1088;
t1138 = t1062 * t1090;
t1134 = pkin(2) * t1159;
t1132 = t1060 * t1160;
t1131 = t1061 * t1160;
t1130 = t1062 * t1160;
t1129 = t1057 * t1160;
t1128 = t1058 * t1160;
t1127 = t1059 * t1160;
t1126 = pkin(3) * t1137;
t1125 = pkin(3) * t1136;
t1124 = pkin(3) * t1135;
t1028 = 0.1e1 / (t1048 * t1091 - t1126);
t1123 = t1028 * t1149;
t1029 = 0.1e1 / (t1048 * t1093 - t1125);
t1122 = t1029 * t1148;
t1030 = 0.1e1 / (t1048 * t1095 - t1124);
t1121 = t1030 * t1147;
t1120 = t1028 * t1158;
t1119 = t1029 * t1157;
t1118 = t1030 * t1156;
t1099 = pkin(3) ^ 2;
t1100 = pkin(2) ^ 2;
t1101 = 0.2e1 * t1073 * t1099 - t1099 + t1100 + 0.2e1 * t1134;
t1049 = pkin(1) * t1160;
t1045 = pkin(1) * t1089 - t1160;
t1044 = pkin(1) * t1087 - t1160;
t1043 = pkin(1) * t1085 - t1160;
t1031 = t1134 + t1100 / 0.2e1 + (t1073 - 0.1e1 / 0.2e1) * t1099;
t1024 = -0.2e1 * t1090 * t1124 + t1111;
t1023 = -0.2e1 * t1088 * t1125 + t1112;
t1022 = -0.2e1 * t1086 * t1126 + t1113;
t1021 = t1101 * t1089 + t1049;
t1020 = t1101 * t1087 + t1049;
t1019 = t1101 * t1085 + t1049;
t1 = [((t1048 * t1138 + t1127) * t1095 + t1089 * (-t1090 * t1130 + t1142)) * t1121 + t1059 * t1153 + ((t1048 * t1139 + t1128) * t1093 + t1087 * (-t1088 * t1131 + t1144)) * t1122 + t1058 * t1154 + ((t1048 * t1140 + t1129) * t1091 + t1085 * (-t1086 * t1132 + t1146)) * t1123 + t1057 * t1155 - g(1) * m(4) + (-((t1031 * t1138 + t1048 * t1127) * t1164 + (t1021 * t1059 + t1024 * t1150) * t1095 - t1062 * t1161 + t1045 * t1142) * t1118 - ((t1031 * t1139 + t1048 * t1128) * t1165 + (t1020 * t1058 + t1023 * t1151) * t1093 - t1061 * t1162 + t1044 * t1144) * t1119 - ((t1031 * t1140 + t1048 * t1129) * t1166 + (t1019 * t1057 + t1022 * t1152) * t1091 - t1060 * t1163 + t1043 * t1146) * t1120) * m(3); ((-t1048 * t1141 + t1130) * t1095 + (t1090 * t1127 + t1150) * t1089) * t1121 + t1062 * t1153 + ((-t1048 * t1143 + t1131) * t1093 + (t1088 * t1128 + t1151) * t1087) * t1122 + t1061 * t1154 + ((-t1048 * t1145 + t1132) * t1091 + (t1086 * t1129 + t1152) * t1085) * t1123 + t1060 * t1155 - g(2) * m(4) + (-((-t1031 * t1141 + t1048 * t1130) * t1164 + (t1062 * t1021 - t1024 * t1142) * t1095 + t1059 * t1161 + t1045 * t1150) * t1118 - ((-t1031 * t1143 + t1048 * t1131) * t1165 + (t1061 * t1020 - t1023 * t1144) * t1093 + t1058 * t1162 + t1044 * t1151) * t1119 - ((-t1031 * t1145 + t1048 * t1132) * t1166 + (t1060 * t1019 - t1022 * t1146) * t1091 + t1057 * t1163 + t1043 * t1152) * t1120) * m(3); -g(3) * m(4) + t1092 * t1149 + t1094 * t1148 + t1096 * t1147 + (-(t1090 * t1072 + (pkin(1) + t1115) * t1096) * t1156 - (t1088 * t1071 + (pkin(1) + t1116) * t1094) * t1157 - (t1086 * t1070 + (pkin(1) + t1117) * t1092) * t1158) * m(3);];
taugX  = t1;
