% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRR1A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*(3+1)/2x8]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRR1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1A0_inertia_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1A0_inertia_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:39
% EndTime: 2019-05-03 14:47:40
% DurationCPUTime: 0.50s
% Computational Cost: add. (473->91), mult. (1057->217), div. (216->8), fcn. (1190->14), ass. (0->95)
t181 = sin(qJ(2,3));
t172 = 0.1e1 / t181;
t187 = xP(3);
t170 = sin(t187);
t171 = cos(t187);
t188 = koppelP(3,2);
t191 = koppelP(3,1);
t157 = t170 * t191 + t171 * t188;
t178 = legFrame(3,3);
t164 = sin(t178);
t167 = cos(t178);
t198 = t170 * t188 - t171 * t191;
t241 = t157 * t167 + t164 * t198;
t244 = t172 * t241;
t182 = sin(qJ(2,2));
t174 = 0.1e1 / t182;
t189 = koppelP(2,2);
t192 = koppelP(2,1);
t158 = t170 * t192 + t171 * t189;
t179 = legFrame(2,3);
t165 = sin(t179);
t168 = cos(t179);
t197 = t170 * t189 - t171 * t192;
t240 = t158 * t168 + t165 * t197;
t243 = t174 * t240;
t183 = sin(qJ(2,1));
t176 = 0.1e1 / t183;
t190 = koppelP(1,2);
t193 = koppelP(1,1);
t159 = t170 * t193 + t171 * t190;
t180 = legFrame(1,3);
t166 = sin(t180);
t169 = cos(t180);
t196 = t170 * t190 - t171 * t193;
t239 = t159 * t169 + t166 * t196;
t242 = t176 * t239;
t194 = 1 / pkin(2);
t238 = 2 * t194;
t184 = cos(qJ(2,3));
t139 = (t164 * t157 - t167 * t198) * t181 - t241 * t184;
t145 = t194 * t244;
t237 = t139 * t145;
t185 = cos(qJ(2,2));
t140 = (t165 * t158 - t168 * t197) * t182 - t240 * t185;
t146 = t194 * t243;
t236 = t140 * t146;
t186 = cos(qJ(2,1));
t141 = (t166 * t159 - t169 * t196) * t183 - t239 * t186;
t147 = t194 * t242;
t235 = t141 * t147;
t151 = -t164 * t181 + t167 * t184;
t152 = t164 * t184 + t181 * t167;
t142 = (-t151 * t157 - t152 * t198) * t172;
t234 = t142 * t241;
t153 = -t165 * t182 + t168 * t185;
t154 = t165 * t185 + t182 * t168;
t143 = (-t153 * t158 - t154 * t197) * t174;
t233 = t143 * t240;
t155 = -t166 * t183 + t169 * t186;
t156 = t166 * t186 + t183 * t169;
t144 = (-t155 * t159 - t156 * t196) * t176;
t232 = t144 * t239;
t231 = t151 * t167;
t230 = t152 * t164;
t173 = 0.1e1 / t181 ^ 2;
t229 = t152 * t173;
t228 = t153 * t168;
t227 = t154 * t165;
t175 = 0.1e1 / t182 ^ 2;
t226 = t154 * t175;
t225 = t155 * t169;
t224 = t156 * t166;
t177 = 0.1e1 / t183 ^ 2;
t223 = t156 * t177;
t219 = t167 * t173;
t218 = t168 * t175;
t217 = t169 * t177;
t216 = t172 * t184;
t215 = t173 * t184;
t214 = t174 * t185;
t213 = t175 * t185;
t212 = t176 * t186;
t211 = t177 * t186;
t210 = t139 * t167 - t151 * t241;
t209 = t139 * t164 - t152 * t241;
t208 = t140 * t168 - t153 * t240;
t207 = t140 * t165 - t154 * t240;
t206 = t141 * t169 - t155 * t239;
t205 = t141 * t166 - t156 * t239;
t204 = t151 * t164 + t152 * t167;
t203 = t153 * t165 + t154 * t168;
t202 = t155 * t166 + t156 * t169;
t195 = 0.1e1 / pkin(2) ^ 2;
t163 = t170 ^ 2 + t171 ^ 2;
t1 = [t151 ^ 2 * t173 + t153 ^ 2 * t175 + t155 ^ 2 * t177, (t167 ^ 2 * t173 + t168 ^ 2 * t175 + t169 ^ 2 * t177) * t195, (-t211 * t225 - t213 * t228 - t215 * t231) * t238, (t172 * t231 + t174 * t228 + t176 * t225) * t238, 0, 0, 0, t163; t151 * t229 + t153 * t226 + t155 * t223, (t164 * t219 + t165 * t218 + t166 * t217) * t195, (-t202 * t211 - t203 * t213 - t204 * t215) * t194, (t172 * t204 + t174 * t203 + t176 * t202) * t194, 0, 0, 0, 0; t152 ^ 2 * t173 + t154 ^ 2 * t175 + t156 ^ 2 * t177, (t164 ^ 2 * t173 + t165 ^ 2 * t175 + t166 ^ 2 * t177) * t195, (-t211 * t224 - t213 * t227 - t215 * t230) * t238, (t172 * t230 + t174 * t227 + t176 * t224) * t238, 0, 0, 0, t163; t139 * t173 * t151 + t140 * t175 * t153 + t141 * t177 * t155, (-t217 * t239 - t218 * t240 - t219 * t241) * t195, (-t206 * t211 - t208 * t213 - t210 * t215) * t194, (t172 * t210 + t174 * t208 + t176 * t206) * t194, 0, -t170, -t171, 0; t139 * t229 + t140 * t226 + t141 * t223, (-t164 * t173 * t241 - t165 * t175 * t240 - t166 * t177 * t239) * t195, (-t205 * t211 - t207 * t213 - t209 * t215) * t194, (t172 * t209 + t174 * t207 + t176 * t205) * t194, 0, t171, -t170, 0; t139 * t172 * t142 + t140 * t174 * t143 + t141 * t176 * t144, (t145 * t244 + t146 * t243 + t147 * t242) * t194, t216 * t237 + t214 * t236 + t212 * t235 + (t212 * t232 + t214 * t233 + t216 * t234) * t194, -t237 - t236 - t235 + (-t232 - t233 - t234) * t194, 1, 0, 0, 0;];
tau_reg  = t1;
