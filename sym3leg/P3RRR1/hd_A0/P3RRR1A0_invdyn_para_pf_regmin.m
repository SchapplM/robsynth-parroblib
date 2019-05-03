% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x10]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RRR1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRR1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRR1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1A0_invdyn_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRR1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1A0_invdyn_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:29
% EndTime: 2019-05-03 15:38:32
% DurationCPUTime: 3.36s
% Computational Cost: add. (25071->268), mult. (39493->514), div. (2880->9), fcn. (29778->26), ass. (0->220)
t205 = 0.1e1 / pkin(1);
t174 = qJ(1,1) + qJ(2,1);
t157 = sin(t174);
t160 = cos(t174);
t186 = sin(qJ(1,1));
t192 = cos(qJ(1,1));
t125 = t157 * t192 - t186 * t160;
t287 = 0.1e1 / t125;
t297 = t287 * t205;
t173 = qJ(1,2) + qJ(2,2);
t156 = sin(t173);
t159 = cos(t173);
t184 = sin(qJ(1,2));
t190 = cos(qJ(1,2));
t124 = t156 * t190 - t184 * t159;
t288 = 0.1e1 / t124;
t296 = t288 * t205;
t172 = qJ(1,3) + qJ(2,3);
t155 = sin(t172);
t158 = cos(t172);
t182 = sin(qJ(1,3));
t188 = cos(qJ(1,3));
t123 = t155 * t188 - t182 * t158;
t289 = 0.1e1 / t123;
t295 = t289 * t205;
t175 = legFrame(3,3);
t161 = sin(t175);
t164 = cos(t175);
t196 = xP(3);
t169 = sin(t196);
t170 = cos(t196);
t197 = koppelP(3,2);
t200 = koppelP(3,1);
t128 = t169 * t200 + t170 * t197;
t193 = xDP(3);
t195 = xDP(1);
t85 = t128 * t193 - t195;
t131 = -t169 * t197 + t170 * t200;
t194 = xDP(2);
t88 = t131 * t193 + t194;
t46 = (t161 * t88 - t85 * t164) * t158 + (t161 * t85 + t88 * t164) * t155;
t176 = legFrame(2,3);
t162 = sin(t176);
t165 = cos(t176);
t198 = koppelP(2,2);
t201 = koppelP(2,1);
t129 = t169 * t201 + t170 * t198;
t86 = t129 * t193 - t195;
t132 = -t169 * t198 + t170 * t201;
t89 = t132 * t193 + t194;
t47 = (t162 * t89 - t86 * t165) * t159 + (t162 * t86 + t89 * t165) * t156;
t177 = legFrame(1,3);
t163 = sin(t177);
t166 = cos(t177);
t199 = koppelP(1,2);
t202 = koppelP(1,1);
t130 = t169 * t202 + t170 * t199;
t87 = t130 * t193 - t195;
t133 = -t169 * t199 + t170 * t202;
t90 = t133 * t193 + t194;
t48 = (t163 * t90 - t87 * t166) * t160 + (t163 * t87 + t90 * t166) * t157;
t140 = -t200 * t182 + t188 * t197;
t143 = t182 * t197 + t188 * t200;
t294 = (t140 * t170 + t169 * t143) * t164 - t161 * (-t169 * t140 + t143 * t170);
t141 = -t201 * t184 + t190 * t198;
t144 = t184 * t198 + t190 * t201;
t293 = (t141 * t170 + t169 * t144) * t165 - t162 * (-t169 * t141 + t144 * t170);
t142 = -t202 * t186 + t192 * t199;
t145 = t186 * t199 + t192 * t202;
t292 = (t142 * t170 + t169 * t145) * t166 - t163 * (-t169 * t142 + t145 * t170);
t291 = -2 * pkin(1);
t290 = 2 * pkin(1);
t203 = pkin(2) ^ 2;
t212 = t155 * t182 + t158 * t188;
t234 = 2 * pkin(2);
t276 = t289 * t46;
t204 = 1 / pkin(2);
t241 = t204 * t205;
t34 = pkin(1) * ((-t182 * t194 - t188 * t195) * t164 - t161 * (-t182 * t195 + t188 * t194) + t294 * t193) - t46 * pkin(2);
t277 = t289 * t34;
t227 = t241 * t277;
t40 = t46 * t295;
t28 = t40 + t227 / 0.2e1;
t31 = t40 + t227;
t286 = (t31 * t203 + (t212 * t28 * t234 + t276) * pkin(1)) * t46;
t211 = t156 * t184 + t159 * t190;
t273 = t288 * t47;
t35 = pkin(1) * ((-t184 * t194 - t190 * t195) * t165 - t162 * (-t184 * t195 + t190 * t194) + t293 * t193) - t47 * pkin(2);
t274 = t288 * t35;
t226 = t241 * t274;
t41 = t47 * t296;
t29 = t41 + t226 / 0.2e1;
t32 = t41 + t226;
t285 = (t32 * t203 + (t211 * t29 * t234 + t273) * pkin(1)) * t47;
t210 = t157 * t186 + t160 * t192;
t270 = t287 * t48;
t36 = pkin(1) * ((-t186 * t194 - t192 * t195) * t166 - t163 * (-t186 * t195 + t192 * t194) + t292 * t193) - t48 * pkin(2);
t271 = t287 * t36;
t225 = t241 * t271;
t42 = t48 * t297;
t30 = t42 + t225 / 0.2e1;
t33 = t42 + t225;
t284 = (t33 * t203 + (t210 * t30 * t234 + t270) * pkin(1)) * t48;
t91 = t164 * t155 + t161 * t158;
t61 = pkin(1) * (t161 * t188 + t164 * t182) + t91 * pkin(2);
t171 = t193 ^ 2;
t178 = xDDP(3);
t179 = xDDP(2);
t73 = -t171 * t128 + t131 * t178 + t179;
t283 = t61 * t73;
t93 = t165 * t156 + t162 * t159;
t62 = pkin(1) * (t162 * t190 + t165 * t184) + t93 * pkin(2);
t74 = -t171 * t129 + t132 * t178 + t179;
t282 = t62 * t74;
t95 = t166 * t157 + t163 * t160;
t63 = pkin(1) * (t163 * t192 + t166 * t186) + t95 * pkin(2);
t75 = -t171 * t130 + t133 * t178 + t179;
t281 = t63 * t75;
t92 = -t155 * t161 + t164 * t158;
t64 = pkin(1) * (-t182 * t161 + t164 * t188) + t92 * pkin(2);
t180 = xDDP(1);
t76 = -t128 * t178 - t171 * t131 + t180;
t280 = t64 * t76;
t94 = -t156 * t162 + t165 * t159;
t65 = pkin(1) * (-t184 * t162 + t165 * t190) + t94 * pkin(2);
t77 = -t129 * t178 - t171 * t132 + t180;
t279 = t65 * t77;
t96 = -t157 * t163 + t166 * t160;
t66 = pkin(1) * (-t186 * t163 + t166 * t192) + t96 * pkin(2);
t78 = -t130 * t178 - t171 * t133 + t180;
t278 = t66 * t78;
t275 = t289 * t91;
t272 = t288 * t93;
t269 = t287 * t95;
t221 = (t161 * t128 + t131 * t164) * t155 - (t128 * t164 - t161 * t131) * t158;
t268 = t289 * (-pkin(1) * t294 + t221 * pkin(2));
t267 = t289 * t221;
t266 = t289 * t61;
t265 = t289 * t64;
t264 = t289 * t92;
t220 = (t162 * t129 + t132 * t165) * t156 - (t129 * t165 - t162 * t132) * t159;
t263 = t288 * (-pkin(1) * t293 + t220 * pkin(2));
t262 = t288 * t220;
t261 = t288 * t62;
t260 = t288 * t65;
t259 = t288 * t94;
t219 = (t163 * t130 + t133 * t166) * t157 - (t130 * t166 - t163 * t133) * t160;
t258 = t287 * (-pkin(1) * t292 + t219 * pkin(2));
t257 = t287 * t219;
t256 = t287 * t63;
t255 = t287 * t66;
t254 = t287 * t96;
t115 = 0.1e1 / t123 ^ 2;
t206 = 1 / pkin(1) ^ 2;
t249 = t115 * t206;
t117 = 0.1e1 / t124 ^ 2;
t247 = t117 * t206;
t119 = 0.1e1 / t125 ^ 2;
t245 = t119 * t206;
t134 = -t161 * g(1) + t164 * g(2);
t137 = t164 * g(1) + t161 * g(2);
t240 = t134 * t155 + t137 * t158;
t135 = -t162 * g(1) + t165 * g(2);
t138 = t165 * g(1) + t162 * g(2);
t239 = t135 * t156 + t138 * t159;
t136 = -t163 * g(1) + t166 * g(2);
t139 = t166 * g(1) + t163 * g(2);
t238 = t136 * t157 + t139 * t160;
t237 = -t134 * t158 + t137 * t155;
t236 = -t135 * t159 + t138 * t156;
t235 = -t136 * t160 + t139 * t157;
t233 = t31 * (t212 * pkin(1) + pkin(2)) * t34;
t232 = t32 * (t211 * pkin(1) + pkin(2)) * t35;
t231 = t33 * (t210 * pkin(1) + pkin(2)) * t36;
t230 = t46 ^ 2 * t249;
t229 = t47 ^ 2 * t247;
t228 = t48 ^ 2 * t245;
t224 = t28 * t227;
t223 = t29 * t226;
t222 = t30 * t225;
t19 = (-(-pkin(2) * t31 - t212 * t276) * t115 * t46 + t31 * t289 * t277) * t206 + (t91 * t73 + t92 * t76) * t295;
t20 = (-(-pkin(2) * t32 - t211 * t273) * t117 * t47 + t32 * t288 * t274) * t206 + (t93 * t74 + t94 * t77) * t296;
t21 = (-(-t33 * pkin(2) - t210 * t270) * t119 * t48 + t33 * t287 * t271) * t206 + (t95 * t75 + t96 * t78) * t297;
t191 = cos(qJ(2,1));
t189 = cos(qJ(2,2));
t187 = cos(qJ(2,3));
t185 = sin(qJ(2,1));
t183 = sin(qJ(2,2));
t181 = sin(qJ(2,3));
t168 = t180 - g(1);
t167 = t179 - g(2);
t127 = -t169 * t178 - t170 * t171;
t126 = -t169 * t171 + t170 * t178;
t101 = t169 * t167 + t170 * t168;
t100 = t170 * t167 - t169 * t168;
t84 = t136 * t186 + t139 * t192;
t83 = -t136 * t192 + t139 * t186;
t82 = t135 * t184 + t138 * t190;
t81 = -t135 * t190 + t138 * t184;
t80 = t134 * t182 + t137 * t188;
t79 = -t134 * t188 + t137 * t182;
t18 = pkin(1) * (-t21 * t185 + t191 * t228) + t238;
t17 = pkin(1) * (t185 * t228 + t21 * t191) + t235;
t16 = pkin(1) * (-t20 * t183 + t189 * t229) + t239;
t15 = pkin(1) * (t183 * t229 + t20 * t189) + t236;
t14 = pkin(1) * (-t19 * t181 + t187 * t230) + t240;
t13 = pkin(1) * (t181 * t230 + t19 * t187) + t237;
t12 = (-(t278 + t281) * t297 + (-t231 - t284) * t245) * t204 + t21;
t11 = (-(t279 + t282) * t296 + (-t232 - t285) * t247) * t204 + t20;
t10 = (-(t280 + t283) * t295 + (-t233 - t286) * t249) * t204 + t19;
t9 = (-(t278 / 0.2e1 + t281 / 0.2e1) * t297 + (-t284 / 0.2e1 - t231 / 0.2e1) * t245) * t204 + t21;
t8 = (-(t279 / 0.2e1 + t282 / 0.2e1) * t296 + (-t285 / 0.2e1 - t232 / 0.2e1) * t247) * t204 + t20;
t7 = (-(t280 / 0.2e1 + t283 / 0.2e1) * t295 + (-t286 / 0.2e1 - t233 / 0.2e1) * t249) * t204 + t19;
t6 = (t185 * t9 + t191 * t222) * t291 + t238;
t5 = (-t185 * t222 + t9 * t191) * t290 + t235;
t4 = (t183 * t8 + t189 * t223) * t291 + t239;
t3 = (-t183 * t223 + t8 * t189) * t290 + t236;
t2 = (t181 * t7 + t187 * t224) * t291 + t240;
t1 = (-t181 * t224 + t7 * t187) * t290 + t237;
t22 = [(t19 * t264 + t20 * t259 + t21 * t254) * t205, (t83 * t254 + t81 * t259 + t79 * t264) * t205, (t84 * t254 + t82 * t259 + t80 * t264) * t205, (t10 * t264 + t11 * t259 + t12 * t254 + (-t10 * t265 - t11 * t260 - t12 * t255) * t204) * t205, (t1 * t264 + t3 * t259 + t5 * t254 + (-t13 * t265 - t15 * t260 - t17 * t255) * t204) * t205, (t2 * t264 + t4 * t259 + t6 * t254 + (-t14 * t265 - t16 * t260 - t18 * t255) * t204) * t205, 0, t127, -t126, -t169 * t100 + t170 * t101; (t19 * t275 + t20 * t272 + t21 * t269) * t205, (t83 * t269 + t81 * t272 + t79 * t275) * t205, (t84 * t269 + t82 * t272 + t80 * t275) * t205, (t10 * t275 + t11 * t272 + t12 * t269 + (-t10 * t266 - t11 * t261 - t12 * t256) * t204) * t205, (t1 * t275 + t3 * t272 + t5 * t269 + (-t13 * t266 - t15 * t261 - t17 * t256) * t204) * t205, (t2 * t275 + t4 * t272 + t6 * t269 + (-t14 * t266 - t16 * t261 - t18 * t256) * t204) * t205, 0, t126, t127, t170 * t100 + t169 * t101; (t19 * t267 + t20 * t262 + t21 * t257) * t205, (t83 * t257 + t81 * t262 + t79 * t267) * t205, (t84 * t257 + t82 * t262 + t80 * t267) * t205, (t10 * t267 + t11 * t262 + t12 * t257 + (-t10 * t268 - t11 * t263 - t12 * t258) * t204) * t205, (t1 * t267 + t3 * t262 + t5 * t257 + (-t13 * t268 - t15 * t263 - t17 * t258) * t204) * t205, (t2 * t267 + t4 * t262 + t6 * t257 + (-t14 * t268 - t16 * t263 - t18 * t258) * t204) * t205, t178, t100, -t101, 0;];
tauX_reg  = t22;
