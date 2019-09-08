% Calculate inertia matrix for parallel robot
% P3PRR1G1P1A0
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
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRR1G1P1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1G1P1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1G1P1A0_inertia_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1G1P1A0_inertia_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PRR1G1P1A0_inertia_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3PRR1G1P1A0_inertia_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'P3PRR1G1P1A0_inertia_para_pf_slag_vp1: Icges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1G1P1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1G1P1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:34
% EndTime: 2019-05-03 14:47:35
% DurationCPUTime: 0.30s
% Computational Cost: add. (517->98), mult. (960->191), div. (324->4), fcn. (1088->14), ass. (0->94)
t181 = sin(qJ(2,3));
t184 = cos(qJ(2,3));
t163 = t184 * rSges(2,1) - t181 * rSges(2,2);
t217 = m(2) * t163;
t182 = sin(qJ(2,2));
t185 = cos(qJ(2,2));
t164 = t185 * rSges(2,1) - t182 * rSges(2,2);
t216 = m(2) * t164;
t183 = sin(qJ(2,1));
t186 = cos(qJ(2,1));
t165 = t186 * rSges(2,1) - t183 * rSges(2,2);
t215 = m(2) * t165;
t197 = 0.1e1 / pkin(2);
t214 = m(2) * t197;
t178 = legFrame(3,3);
t167 = sin(t178);
t170 = cos(t178);
t149 = -t167 * t181 + t170 * t184;
t175 = 0.1e1 / t181;
t213 = t149 * t175;
t150 = t167 * t184 + t181 * t170;
t212 = t150 * t175;
t179 = legFrame(2,3);
t168 = sin(t179);
t171 = cos(t179);
t151 = -t168 * t182 + t171 * t185;
t176 = 0.1e1 / t182;
t211 = t151 * t176;
t152 = t168 * t185 + t182 * t171;
t210 = t152 * t176;
t180 = legFrame(1,3);
t169 = sin(t180);
t172 = cos(t180);
t153 = -t169 * t183 + t172 * t186;
t177 = 0.1e1 / t183;
t209 = t153 * t177;
t154 = t169 * t186 + t183 * t172;
t208 = t154 * t177;
t166 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3);
t207 = t166 * t197;
t206 = t167 * t175;
t205 = t168 * t176;
t204 = t169 * t177;
t203 = t170 * t175;
t202 = t171 * t176;
t201 = t172 * t177;
t200 = t163 * t214;
t199 = t164 * t214;
t198 = t165 * t214;
t196 = koppelP(1,1);
t195 = koppelP(2,1);
t194 = koppelP(3,1);
t193 = koppelP(1,2);
t192 = koppelP(2,2);
t191 = koppelP(3,2);
t190 = rSges(3,1);
t189 = rSges(3,2);
t188 = xP(3);
t187 = m(1) + m(2);
t174 = cos(t188);
t173 = sin(t188);
t162 = -t173 * t193 + t174 * t196;
t161 = -t173 * t192 + t174 * t195;
t160 = -t173 * t191 + t174 * t194;
t159 = -t173 * t196 - t174 * t193;
t158 = -t173 * t195 - t174 * t192;
t157 = -t173 * t194 - t174 * t191;
t156 = m(3) * (-t173 * t189 + t174 * t190);
t155 = m(3) * (-t173 * t190 - t174 * t189);
t148 = (t154 * t187 - t169 * t198) * t177;
t147 = (t153 * t187 - t172 * t198) * t177;
t146 = (t152 * t187 - t168 * t199) * t176;
t145 = (t151 * t187 - t171 * t199) * t176;
t144 = (t150 * t187 - t167 * t200) * t175;
t143 = (t149 * t187 - t170 * t200) * t175;
t142 = (-t159 * t172 - t162 * t169) * t197 * t177;
t141 = (-t158 * t171 - t161 * t168) * t197 * t176;
t140 = (-t157 * t170 - t160 * t167) * t197 * t175;
t139 = (t154 * t215 - t169 * t207) * t177;
t138 = (t153 * t215 - t172 * t207) * t177;
t137 = (t152 * t216 - t168 * t207) * t176;
t136 = (t151 * t216 - t171 * t207) * t176;
t135 = (t150 * t217 - t167 * t207) * t175;
t134 = (t149 * t217 - t170 * t207) * t175;
t133 = (t153 * t159 + t154 * t162) * t177;
t132 = (t151 * t158 + t152 * t161) * t176;
t131 = (t149 * t157 + t150 * t160) * t175;
t130 = t133 * t187 + t142 * t215;
t129 = t132 * t187 + t141 * t216;
t128 = t131 * t187 + t140 * t217;
t127 = t133 * t215 + t142 * t166;
t126 = t132 * t216 + t141 * t166;
t125 = t131 * t217 + t140 * t166;
t1 = [t143 * t213 + t145 * t211 + t147 * t209 + m(3) + (-t134 * t203 - t136 * t202 - t138 * t201) * t197, t143 * t212 + t145 * t210 + t147 * t208 + (-t134 * t206 - t136 * t205 - t138 * t204) * t197, t143 * t131 + t145 * t132 + t147 * t133 + t134 * t140 + t136 * t141 + t138 * t142 + t155; t144 * t213 + t146 * t211 + t148 * t209 + (-t135 * t203 - t137 * t202 - t139 * t201) * t197, t144 * t212 + t146 * t210 + t148 * t208 + m(3) + (-t135 * t206 - t137 * t205 - t139 * t204) * t197, t144 * t131 + t146 * t132 + t148 * t133 + t135 * t140 + t137 * t141 + t139 * t142 + t156; t128 * t213 + t129 * t211 + t130 * t209 + t155 + (-t125 * t203 - t126 * t202 - t127 * t201) * t197, t128 * t212 + t129 * t210 + t130 * t208 + t156 + (-t125 * t206 - t126 * t205 - t127 * t204) * t197, t130 * t133 + t127 * t142 + t129 * t132 + t126 * t141 + t128 * t131 + t125 * t140 + Icges(3,3) + m(3) * (t189 ^ 2 + t190 ^ 2);];
MX  = t1;
