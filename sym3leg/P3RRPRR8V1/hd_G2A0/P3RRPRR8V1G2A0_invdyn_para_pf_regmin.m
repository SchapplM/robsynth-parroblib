% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRPRR8V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
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
%   pkin=[a3,a4,d1,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x13]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RRPRR8V1G2A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:59:24
% EndTime: 2020-08-06 19:59:31
% DurationCPUTime: 7.04s
% Computational Cost: add. (21838->429), mult. (47553->749), div. (2640->21), fcn. (32901->35), ass. (0->344)
t196 = xDP(2);
t171 = sin(pkin(5));
t404 = pkin(2) * t171;
t128 = t196 * t404;
t197 = xDP(1);
t129 = t197 * t404;
t176 = legFrame(3,2);
t142 = sin(t176);
t145 = cos(t176);
t182 = sin(qJ(2,3));
t183 = sin(qJ(1,3));
t189 = cos(qJ(1,3));
t195 = xDP(3);
t313 = t195 * t404;
t172 = cos(pkin(5));
t394 = t172 * pkin(2);
t132 = pkin(1) + t394;
t421 = t197 * t132;
t422 = t196 * t132;
t425 = t182 * ((-t183 * t129 + t422) * t145 + (t183 * t128 + t421) * t142 - t189 * t313);
t177 = legFrame(2,2);
t143 = sin(t177);
t146 = cos(t177);
t184 = sin(qJ(2,2));
t185 = sin(qJ(1,2));
t191 = cos(qJ(1,2));
t424 = t184 * ((-t185 * t129 + t422) * t146 + (t185 * t128 + t421) * t143 - t191 * t313);
t178 = legFrame(1,2);
t144 = sin(t178);
t147 = cos(t178);
t186 = sin(qJ(2,1));
t187 = sin(qJ(1,1));
t193 = cos(qJ(1,1));
t423 = t186 * ((-t187 * t129 + t422) * t147 + (t187 * t128 + t421) * t144 - t193 * t313);
t420 = 0.2e1 * pkin(1);
t173 = pkin(4) + qJ(3,3);
t209 = (t173 ^ 2);
t160 = 0.1e1 / t209;
t133 = t183 * t173;
t188 = cos(qJ(2,3));
t326 = t189 * t173;
t342 = t132 * t188;
t34 = ((t183 * t421 + t128) * t188 - t197 * t326) * t145 + ((-t183 * t422 + t129) * t188 + t196 * t326) * t142 + t195 * (t189 * t342 + t133) + t425;
t339 = t132 * t195;
t40 = (t189 * t339 + t128 * t145 + t129 * t142 + (-t142 * t422 + t145 * t421) * t183) * t188 + t425;
t329 = t171 * t182;
t219 = pkin(2) * t329 - t342;
t88 = 0.1e1 / t219;
t387 = t40 * t88;
t316 = t34 * t387;
t266 = t160 * t316;
t201 = pkin(2) ^ 2;
t202 = pkin(1) ^ 2;
t320 = -t201 - t202;
t124 = t394 * t420 - t320;
t148 = t188 * pkin(1);
t159 = 0.1e1 / t173;
t179 = xDDP(3);
t180 = xDDP(2);
t181 = xDDP(1);
t161 = t159 * t160;
t265 = t161 * t316;
t103 = t142 * t197 + t145 * t196;
t156 = qJ(2,3) + pkin(5);
t416 = t148 + pkin(2) * cos(t156);
t360 = t103 ^ 2 / t416 ^ 2;
t291 = t142 * t404;
t294 = t145 * t404;
t338 = t142 * t183;
t58 = (-t132 * t338 + t294) * t188 + t182 * (t132 * t145 + t183 * t291);
t335 = t145 * t183;
t61 = (t132 * t335 + t291) * t188 + (t142 * t132 - t183 * t294) * t182;
t22 = t265 + (t189 * t179 + t124 / (t148 + (t172 * t188 - t329) * pkin(2)) * t360 - (t58 * t180 + t61 * t181 + (t219 * t387 - t34) * t160 * t40) * t88) * t159;
t382 = qJ(3,3) * t22;
t419 = 0.2e1 * t266 - t382;
t174 = pkin(4) + qJ(3,2);
t211 = (t174 ^ 2);
t163 = 0.1e1 / t211;
t134 = t185 * t174;
t190 = cos(qJ(2,2));
t325 = t191 * t174;
t341 = t132 * t190;
t35 = ((t185 * t421 + t128) * t190 - t197 * t325) * t146 + ((-t185 * t422 + t129) * t190 + t196 * t325) * t143 + t195 * (t191 * t341 + t134) + t424;
t41 = (t191 * t339 + t128 * t146 + t129 * t143 + (-t143 * t422 + t146 * t421) * t185) * t190 + t424;
t328 = t171 * t184;
t218 = pkin(2) * t328 - t341;
t90 = 0.1e1 / t218;
t386 = t41 * t90;
t315 = t35 * t386;
t264 = t163 * t315;
t149 = t190 * pkin(1);
t162 = 0.1e1 / t174;
t164 = t162 * t163;
t263 = t164 * t315;
t104 = t143 * t197 + t146 * t196;
t157 = qJ(2,2) + pkin(5);
t415 = t149 + pkin(2) * cos(t157);
t359 = t104 ^ 2 / t415 ^ 2;
t290 = t143 * t404;
t293 = t146 * t404;
t337 = t143 * t185;
t59 = (-t132 * t337 + t293) * t190 + t184 * (t132 * t146 + t185 * t290);
t334 = t146 * t185;
t62 = (t132 * t334 + t290) * t190 + (t143 * t132 - t185 * t293) * t184;
t23 = t263 + (t191 * t179 + t124 / (t149 + (t172 * t190 - t328) * pkin(2)) * t359 - (t59 * t180 + t62 * t181 + (t218 * t386 - t35) * t163 * t41) * t90) * t162;
t383 = qJ(3,2) * t23;
t418 = 0.2e1 * t264 - t383;
t175 = pkin(4) + qJ(3,1);
t213 = (t175 ^ 2);
t166 = 0.1e1 / t213;
t135 = t187 * t175;
t192 = cos(qJ(2,1));
t324 = t193 * t175;
t340 = t132 * t192;
t36 = ((t187 * t421 + t128) * t192 - t197 * t324) * t147 + ((-t187 * t422 + t129) * t192 + t196 * t324) * t144 + t195 * (t193 * t340 + t135) + t423;
t42 = (t193 * t339 + t128 * t147 + t129 * t144 + (-t144 * t422 + t147 * t421) * t187) * t192 + t423;
t327 = t171 * t186;
t217 = pkin(2) * t327 - t340;
t92 = 0.1e1 / t217;
t385 = t42 * t92;
t314 = t36 * t385;
t262 = t166 * t314;
t150 = t192 * pkin(1);
t165 = 0.1e1 / t175;
t167 = t165 * t166;
t261 = t167 * t314;
t105 = t144 * t197 + t147 * t196;
t158 = qJ(2,1) + pkin(5);
t414 = t150 + pkin(2) * cos(t158);
t358 = t105 ^ 2 / t414 ^ 2;
t289 = t144 * t404;
t292 = t147 * t404;
t336 = t144 * t187;
t60 = (-t132 * t336 + t292) * t192 + t186 * (t132 * t147 + t187 * t289);
t333 = t147 * t187;
t63 = (t132 * t333 + t289) * t192 + (t144 * t132 - t187 * t292) * t186;
t24 = t261 + (t193 * t179 + t124 / (t150 + (t172 * t192 - t327) * pkin(2)) * t358 - (t60 * t180 + t63 * t181 + (t217 * t385 - t36) * t166 * t42) * t92) * t165;
t384 = qJ(3,1) * t24;
t417 = 0.2e1 * t262 - t384;
t397 = t145 * g(1);
t400 = t142 * g(2);
t231 = t397 - t400;
t403 = g(3) * t183;
t82 = -t189 * t231 + t403;
t396 = t146 * g(1);
t399 = t143 * g(2);
t230 = t396 - t399;
t402 = g(3) * t185;
t83 = -t191 * t230 + t402;
t395 = t147 * g(1);
t398 = t144 * g(2);
t229 = t395 - t398;
t401 = g(3) * t187;
t84 = -t193 * t229 + t401;
t115 = 0.1e1 / t416;
t117 = 0.1e1 / t415;
t119 = 0.1e1 / t414;
t168 = t188 ^ 2;
t240 = t416 * t265;
t393 = t182 * pkin(1);
t237 = pkin(2) * sin(t156) + t393;
t357 = t103 * t115;
t25 = -t124 * t357 - t237 * t387;
t198 = 0.2e1 * qJ(2,3);
t276 = pkin(2) * t420;
t269 = -(0.2e1 * t237 * t173 * t357 + (-(-t201 * cos(0.2e1 * t156) - cos(t198) * t202 - (2 * t209) + (-cos(t198 + pkin(5)) - t172) * t276 + t320) * t387 + 0.2e1 * t416 * t34) * t159) * t387 / 0.2e1;
t272 = t387 * t393;
t73 = t189 * t416 + t133;
t363 = t73 * t179;
t70 = -t183 * t219 - t326;
t94 = t182 * t132 + t188 * t404;
t53 = -t70 * t142 + t94 * t145;
t368 = t53 * t180;
t52 = t94 * t142 + t70 * t145;
t369 = t52 * t181;
t49 = (t142 * t181 + t145 * t180 + t237 * t360) * t115;
t372 = t49 * t182;
t151 = g(3) * t189;
t85 = t231 * t183 + t151;
t1 = t168 * t202 * t22 + (-pkin(1) * t372 - t419 - t85) * qJ(3,3) + (-0.2e1 * (t397 / 0.2e1 - t400 / 0.2e1) * t189 + t403 + t160 * t269 - t240 - 0.2e1 * (t369 / 0.2e1 + t368 / 0.2e1 + t363 / 0.2e1 + (-t272 - t25 / 0.2e1) * t357) * t159 - qJ(3,3) * t360) * t148;
t413 = t1 * t88;
t169 = t190 ^ 2;
t239 = t415 * t263;
t392 = t184 * pkin(1);
t236 = pkin(2) * sin(t157) + t392;
t356 = t104 * t117;
t26 = -t124 * t356 - t236 * t386;
t199 = 0.2e1 * qJ(2,2);
t268 = -(0.2e1 * t236 * t174 * t356 + (-(-t201 * cos(0.2e1 * t157) - cos(t199) * t202 - (2 * t211) + (-cos(pkin(5) + t199) - t172) * t276 + t320) * t386 + 0.2e1 * t415 * t35) * t162) * t386 / 0.2e1;
t271 = t386 * t392;
t74 = t191 * t415 + t134;
t362 = t74 * t179;
t71 = -t185 * t218 - t325;
t95 = t184 * t132 + t190 * t404;
t55 = -t71 * t143 + t95 * t146;
t366 = t55 * t180;
t54 = t95 * t143 + t71 * t146;
t367 = t54 * t181;
t50 = (t143 * t181 + t146 * t180 + t236 * t359) * t117;
t371 = t50 * t184;
t152 = g(3) * t191;
t86 = t230 * t185 + t152;
t2 = t169 * t202 * t23 + (-pkin(1) * t371 - t418 - t86) * qJ(3,2) + (-0.2e1 * (t396 / 0.2e1 - t399 / 0.2e1) * t191 + t402 + t163 * t268 - t239 - 0.2e1 * (t367 / 0.2e1 + t366 / 0.2e1 + t362 / 0.2e1 + (-t271 - t26 / 0.2e1) * t356) * t162 - qJ(3,2) * t359) * t149;
t412 = t2 * t90;
t170 = t192 ^ 2;
t238 = t414 * t261;
t200 = 0.2e1 * qJ(2,1);
t391 = t186 * pkin(1);
t235 = pkin(2) * sin(t158) + t391;
t355 = t105 * t119;
t267 = -(0.2e1 * t235 * t175 * t355 + (-(-t201 * cos(0.2e1 * t158) - cos(t200) * t202 - (2 * t213) + (-cos(pkin(5) + t200) - t172) * t276 + t320) * t385 + 0.2e1 * t414 * t36) * t165) * t385 / 0.2e1;
t27 = -t124 * t355 - t235 * t385;
t270 = t385 * t391;
t75 = t193 * t414 + t135;
t361 = t75 * t179;
t72 = -t187 * t217 - t324;
t96 = t186 * t132 + t192 * t404;
t57 = -t72 * t144 + t96 * t147;
t364 = t57 * t180;
t56 = t96 * t144 + t72 * t147;
t365 = t56 * t181;
t51 = (t144 * t181 + t147 * t180 + t235 * t358) * t119;
t370 = t51 * t186;
t153 = g(3) * t193;
t87 = t229 * t187 + t153;
t3 = t170 * t202 * t24 + (-pkin(1) * t370 - t417 - t87) * qJ(3,1) + (-0.2e1 * (t395 / 0.2e1 - t398 / 0.2e1) * t193 + t401 + t166 * t267 - t238 - 0.2e1 * (t365 / 0.2e1 + t364 / 0.2e1 + t361 / 0.2e1 + (-t270 - t27 / 0.2e1) * t355) * t165 - qJ(3,1) * t358) * t150;
t411 = t3 * t92;
t410 = 0.2e1 * t168 - 0.1e1;
t409 = 0.2e1 * t169 - 0.1e1;
t408 = 0.2e1 * t170 - 0.1e1;
t390 = t40 ^ 2 / t219 ^ 2;
t389 = t41 ^ 2 / t218 ^ 2;
t388 = t42 ^ 2 / t217 ^ 2;
t381 = t159 * t88;
t380 = t162 * t90;
t379 = t165 * t92;
t378 = t188 * t22;
t377 = t190 * t23;
t376 = t192 * t24;
t375 = t22 * t182;
t374 = t23 * t184;
t373 = t24 * t186;
t112 = t142 * g(1) + t145 * g(2);
t354 = t112 * t188;
t113 = t143 * g(1) + t146 * g(2);
t353 = t113 * t190;
t114 = t144 * g(1) + t147 * g(2);
t352 = t114 * t192;
t351 = t115 * t142;
t350 = t115 * t145;
t349 = t115 * t182;
t348 = t117 * t143;
t347 = t117 * t146;
t346 = t117 * t184;
t345 = t119 * t144;
t344 = t119 * t147;
t343 = t119 * t186;
t332 = t159 * t189;
t331 = t162 * t191;
t330 = t165 * t193;
t312 = t58 * t381;
t311 = t61 * t381;
t310 = t82 * t381;
t309 = t85 * t381;
t308 = t160 * t390;
t307 = t161 * t390;
t306 = t59 * t380;
t305 = t62 * t380;
t304 = t83 * t380;
t303 = t86 * t380;
t302 = t163 * t389;
t301 = t164 * t389;
t300 = t60 * t379;
t299 = t63 * t379;
t298 = t84 * t379;
t297 = t87 * t379;
t296 = t166 * t388;
t295 = t167 * t388;
t288 = t22 * t349;
t287 = t115 * t378;
t286 = t23 * t346;
t285 = t117 * t377;
t284 = t24 * t343;
t283 = t119 * t376;
t282 = t159 * t188 * t82;
t281 = t162 * t190 * t83;
t280 = t165 * t192 * t84;
t279 = t82 * t332;
t278 = t83 * t331;
t277 = t84 * t330;
t228 = t40 * t357 * t381;
t13 = t188 * t375 - t410 * t228;
t275 = -0.2e1 * t13 * t381;
t227 = t41 * t356 * t380;
t14 = t190 * t374 - t409 * t227;
t274 = -0.2e1 * t14 * t380;
t226 = t42 * t355 * t379;
t15 = t192 * t373 - t408 * t226;
t273 = -0.2e1 * t15 * t379;
t254 = t182 * t310;
t253 = t88 * t282;
t252 = t188 * t308;
t251 = t184 * t304;
t250 = t90 * t281;
t249 = t190 * t302;
t248 = t186 * t298;
t247 = t92 * t280;
t246 = t192 * t296;
t225 = -g(1) * t335 + g(2) * t338 - t151;
t224 = -g(1) * t334 + g(2) * t337 - t152;
t223 = -g(1) * t333 + g(2) * t336 - t153;
t222 = t252 * t349;
t221 = t249 * t346;
t220 = t246 * t343;
t43 = t188 * t360 + t372;
t45 = t190 * t359 + t371;
t47 = t192 * t358 + t370;
t216 = -t142 * t288 - t143 * t286 - t144 * t284;
t215 = -t145 * t288 - t146 * t286 - t147 * t284;
t69 = t114 * t186 - t192 * t223;
t68 = -t186 * t223 - t352;
t67 = t113 * t184 - t190 * t224;
t66 = -t184 * t224 - t353;
t65 = t112 * t182 - t188 * t225;
t64 = -t182 * t225 - t354;
t48 = -t186 * t358 + t51 * t192;
t46 = -t184 * t359 + t50 * t190;
t44 = -t182 * t360 + t49 * t188;
t30 = t408 * t296;
t29 = t409 * t302;
t28 = t410 * t308;
t18 = (-0.2e1 * t192 * t226 + t373) * t186;
t17 = (-0.2e1 * t190 * t227 + t374) * t184;
t16 = (-0.2e1 * t188 * t228 + t375) * t182;
t12 = -pkin(1) * t47 + t223 - 0.2e1 * t262 + 0.2e1 * t384;
t11 = -pkin(1) * t45 + t224 - 0.2e1 * t264 + 0.2e1 * t383;
t10 = -pkin(1) * t43 + t225 - 0.2e1 * t266 + 0.2e1 * t382;
t9 = (pkin(1) * t246 - t223 + t417) * t186 - t352 + t51 * pkin(1);
t8 = (pkin(1) * t249 - t224 + t418) * t184 - t353 + t50 * pkin(1);
t7 = (pkin(1) * t252 - t225 + t419) * t182 - t354 + t49 * pkin(1);
t6 = t238 - pkin(1) * t376 + (-qJ(3,1) * t388 - t267) * t166 + (t361 + t364 + t365 + (-t27 - 0.2e1 * t270) * t355) * t165 - t84;
t5 = t239 - pkin(1) * t377 + (-qJ(3,2) * t389 - t268) * t163 + (t362 + t366 + t367 + (-t26 - 0.2e1 * t271) * t356) * t162 - t83;
t4 = t240 - pkin(1) * t378 + (-qJ(3,3) * t390 - t269) * t160 + (t363 + t368 + t369 + (-t25 - 0.2e1 * t272) * t357) * t159 - t82;
t19 = [-t22 * t311 - t23 * t305 - t24 * t299, -t298 * t63 - t304 * t62 - t310 * t61, -t297 * t63 - t303 * t62 - t309 * t61, -t142 * t222 - t143 * t221 - t144 * t220 - t16 * t311 - t17 * t305 - t18 * t299, t273 * t63 + t274 * t62 + t275 * t61 - t28 * t351 - t29 * t348 - t30 * t345, -t299 * t47 - t305 * t45 - t311 * t43 - t216, t142 * t287 + t143 * t285 + t144 * t283 - t299 * t48 - t305 * t46 - t311 * t44, t51 * t345 + t50 * t348 + t49 * t351, -t247 * t63 - t250 * t62 - t253 * t61 + t68 * t345 + t66 * t348 + t64 * t351, t248 * t63 + t251 * t62 + t254 * t61 + t345 * t69 + t348 * t67 + t351 * t65, pkin(1) * t216 - t10 * t311 - t11 * t305 - t12 * t299 - t295 * t56 - t301 * t54 - t307 * t52, (-t63 * t411 + t56 * t6) * t165 + (-t62 * t412 + t5 * t54) * t162 + (t4 * t52 - t61 * t413) * t159 + (t345 * t9 + t348 * t8 + t351 * t7) * pkin(1), t181 - g(1); -t22 * t312 - t23 * t306 - t24 * t300, -t298 * t60 - t304 * t59 - t310 * t58, -t297 * t60 - t303 * t59 - t309 * t58, -t145 * t222 - t146 * t221 - t147 * t220 - t16 * t312 - t17 * t306 - t18 * t300, t273 * t60 + t274 * t59 + t275 * t58 - t28 * t350 - t29 * t347 - t30 * t344, -t300 * t47 - t306 * t45 - t312 * t43 - t215, t145 * t287 + t146 * t285 + t147 * t283 - t300 * t48 - t306 * t46 - t312 * t44, t51 * t344 + t50 * t347 + t49 * t350, -t247 * t60 - t250 * t59 - t253 * t58 + t68 * t344 + t66 * t347 + t64 * t350, t248 * t60 + t251 * t59 + t254 * t58 + t344 * t69 + t347 * t67 + t350 * t65, pkin(1) * t215 - t10 * t312 - t11 * t306 - t12 * t300 - t295 * t57 - t301 * t55 - t307 * t53, (-t60 * t411 + t57 * t6) * t165 + (-t59 * t412 + t5 * t55) * t162 + (t4 * t53 - t58 * t413) * t159 + (t344 * t9 + t347 * t8 + t350 * t7) * pkin(1), t180 - g(2); t22 * t332 + t23 * t331 + t24 * t330, t277 + t278 + t279, t87 * t330 + t86 * t331 + t85 * t332, t16 * t332 + t17 * t331 + t18 * t330, 0.2e1 * t13 * t332 + 0.2e1 * t14 * t331 + 0.2e1 * t15 * t330, t47 * t330 + t45 * t331 + t43 * t332, t48 * t330 + t46 * t331 + t44 * t332, 0, t189 * t282 + t191 * t281 + t193 * t280, -t182 * t279 - t184 * t278 - t186 * t277, t10 * t332 + t11 * t331 + t12 * t330 - t295 * t75 - t301 * t74 - t307 * t73, (t193 * t3 + t6 * t75) * t165 + (t191 * t2 + t5 * t74) * t162 + (t1 * t189 + t4 * t73) * t159, t179 - g(3);];
tauX_reg  = t19;
