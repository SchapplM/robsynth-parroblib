% Calculate inertia matrix for parallel robot
% P3PRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRR1G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1G1A0_inertia_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1G1A0_inertia_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PRR1G1A0_inertia_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3PRR1G1A0_inertia_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'P3PRR1G1A0_inertia_para_pf_slag_vp2: Ifges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:37
% EndTime: 2019-05-03 14:47:37
% DurationCPUTime: 0.25s
% Computational Cost: add. (461->94), mult. (818->181), div. (324->4), fcn. (1088->14), ass. (0->91)
t168 = legFrame(3,3);
t157 = sin(t168);
t160 = cos(t168);
t171 = sin(qJ(2,3));
t174 = cos(qJ(2,3));
t141 = -t157 * t171 + t160 * t174;
t165 = 0.1e1 / t171;
t207 = t141 * t165;
t142 = t157 * t174 + t171 * t160;
t206 = t142 * t165;
t169 = legFrame(2,3);
t158 = sin(t169);
t161 = cos(t169);
t172 = sin(qJ(2,2));
t175 = cos(qJ(2,2));
t143 = -t158 * t172 + t161 * t175;
t166 = 0.1e1 / t172;
t205 = t143 * t166;
t144 = t158 * t175 + t172 * t161;
t204 = t144 * t166;
t170 = legFrame(1,3);
t159 = sin(t170);
t162 = cos(t170);
t173 = sin(qJ(2,1));
t176 = cos(qJ(2,1));
t145 = -t159 * t173 + t162 * t176;
t167 = 0.1e1 / t173;
t203 = t145 * t167;
t146 = t159 * t176 + t173 * t162;
t202 = t146 * t167;
t201 = t157 * t165;
t187 = 0.1e1 / pkin(2);
t200 = t157 * t187;
t199 = t158 * t166;
t198 = t158 * t187;
t197 = t159 * t167;
t196 = t159 * t187;
t195 = t160 * t165;
t194 = t160 * t187;
t193 = t161 * t166;
t192 = t161 * t187;
t191 = t162 * t167;
t190 = t162 * t187;
t178 = xP(3);
t163 = sin(t178);
t164 = cos(t178);
t179 = mrSges(3,2);
t180 = mrSges(3,1);
t189 = -t163 * t179 + t164 * t180;
t188 = -t163 * t180 - t164 * t179;
t186 = koppelP(1,1);
t185 = koppelP(2,1);
t184 = koppelP(3,1);
t183 = koppelP(1,2);
t182 = koppelP(2,2);
t181 = koppelP(3,2);
t177 = m(1) + m(2);
t155 = t176 * mrSges(2,1) - t173 * mrSges(2,2);
t154 = t175 * mrSges(2,1) - t172 * mrSges(2,2);
t153 = t174 * mrSges(2,1) - t171 * mrSges(2,2);
t152 = -t163 * t183 + t164 * t186;
t151 = -t163 * t182 + t164 * t185;
t150 = -t163 * t181 + t164 * t184;
t149 = -t163 * t186 - t164 * t183;
t148 = -t163 * t185 - t164 * t182;
t147 = -t163 * t184 - t164 * t181;
t140 = (-Ifges(2,3) * t196 + t146 * t155) * t167;
t139 = (-Ifges(2,3) * t190 + t145 * t155) * t167;
t138 = (-Ifges(2,3) * t198 + t144 * t154) * t166;
t137 = (-Ifges(2,3) * t192 + t143 * t154) * t166;
t136 = (-Ifges(2,3) * t200 + t142 * t153) * t165;
t135 = (-Ifges(2,3) * t194 + t141 * t153) * t165;
t134 = (t146 * t177 - t155 * t196) * t167;
t133 = (t145 * t177 - t155 * t190) * t167;
t132 = (t144 * t177 - t154 * t198) * t166;
t131 = (t143 * t177 - t154 * t192) * t166;
t130 = (t142 * t177 - t153 * t200) * t165;
t129 = (t141 * t177 - t153 * t194) * t165;
t128 = (-t149 * t162 - t152 * t159) * t187 * t167;
t127 = (-t148 * t161 - t151 * t158) * t187 * t166;
t126 = (-t147 * t160 - t150 * t157) * t187 * t165;
t125 = (t145 * t149 + t146 * t152) * t167;
t124 = (t143 * t148 + t144 * t151) * t166;
t123 = (t141 * t147 + t142 * t150) * t165;
t122 = t128 * Ifges(2,3) + t125 * t155;
t121 = t127 * Ifges(2,3) + t124 * t154;
t120 = t126 * Ifges(2,3) + t123 * t153;
t119 = t125 * t177 + t128 * t155;
t118 = t124 * t177 + t127 * t154;
t117 = t123 * t177 + t126 * t153;
t1 = [t129 * t207 + t131 * t205 + t133 * t203 + m(3) + (-t135 * t195 - t137 * t193 - t139 * t191) * t187, t129 * t206 + t131 * t204 + t133 * t202 + (-t135 * t201 - t137 * t199 - t139 * t197) * t187, t129 * t123 + t131 * t124 + t133 * t125 + t135 * t126 + t137 * t127 + t139 * t128 + t188; t130 * t207 + t132 * t205 + t134 * t203 + (-t136 * t195 - t138 * t193 - t140 * t191) * t187, t130 * t206 + t132 * t204 + t134 * t202 + m(3) + (-t136 * t201 - t138 * t199 - t140 * t197) * t187, t130 * t123 + t132 * t124 + t134 * t125 + t136 * t126 + t138 * t127 + t140 * t128 + t189; t117 * t207 + t118 * t205 + t119 * t203 + (-t120 * t195 - t121 * t193 - t122 * t191) * t187 + t188, t117 * t206 + t118 * t204 + t119 * t202 + (-t120 * t201 - t121 * t199 - t122 * t197) * t187 + t189, t117 * t123 + t118 * t124 + t119 * t125 + t120 * t126 + t121 * t127 + t122 * t128 + Ifges(3,3);];
MX  = t1;
