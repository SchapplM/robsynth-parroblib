% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPRRR6V1G3A0
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
% tauX_reg [3x12]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPRRR6V1G3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:42:01
% EndTime: 2020-08-06 18:42:14
% DurationCPUTime: 13.07s
% Computational Cost: add. (99660->553), mult. (115080->943), div. (16212->22), fcn. (82230->84), ass. (0->408)
t239 = sin(pkin(7));
t265 = -pkin(6) - pkin(5);
t169 = t265 * t239;
t156 = -t169 + pkin(1);
t255 = cos(qJ(1,3));
t125 = t156 * t255;
t240 = cos(pkin(7));
t249 = sin(qJ(1,3));
t403 = t249 * t265;
t404 = t249 * t239;
t508 = t125 - pkin(2) * t404 + (pkin(2) * t255 - t403) * t240;
t257 = cos(qJ(1,2));
t126 = t156 * t257;
t251 = sin(qJ(1,2));
t400 = t251 * t265;
t401 = t251 * t239;
t507 = t126 - pkin(2) * t401 + (pkin(2) * t257 - t400) * t240;
t259 = cos(qJ(1,1));
t127 = t156 * t259;
t253 = sin(qJ(1,1));
t397 = t253 * t265;
t398 = t253 * t239;
t506 = t127 - pkin(2) * t398 + (pkin(2) * t259 - t397) * t240;
t242 = legFrame(3,2);
t191 = sin(t242);
t194 = cos(t242);
t122 = t194 * g(1) - t191 * g(2);
t91 = g(3) * t255 + t122 * t249;
t92 = -g(3) * t249 + t122 * t255;
t243 = legFrame(2,2);
t192 = sin(t243);
t195 = cos(t243);
t123 = t195 * g(1) - t192 * g(2);
t93 = g(3) * t257 + t123 * t251;
t94 = -g(3) * t251 + t123 * t257;
t244 = legFrame(1,2);
t193 = sin(t244);
t196 = cos(t244);
t124 = t196 * g(1) - t193 * g(2);
t95 = g(3) * t259 + t124 * t253;
t96 = -g(3) * t253 + t124 * t259;
t283 = 0.2e1 * pkin(2);
t213 = qJ(1,1) + pkin(7);
t164 = t244 + t213;
t165 = -t244 + t213;
t108 = -sin(t164) - sin(t165);
t111 = cos(t165) - cos(t164);
t258 = cos(qJ(3,1));
t199 = t258 * pkin(3);
t178 = t199 + pkin(2);
t190 = t240 * pkin(1);
t142 = t190 + t178;
t134 = 0.1e1 / t142;
t135 = 0.1e1 / t142 ^ 2;
t189 = cos(t213);
t246 = xDDP(2);
t247 = xDDP(1);
t252 = sin(qJ(3,1));
t456 = pkin(3) * t252;
t100 = t240 * t178 + t156;
t168 = t265 * t240;
t406 = t239 * t178;
t104 = t168 + t406;
t262 = xDP(3);
t263 = xDP(2);
t264 = xDP(1);
t286 = 0.1e1 / pkin(3);
t166 = t190 + pkin(2);
t399 = t252 * t166;
t234 = 0.1e1 / t252;
t413 = t134 * t234;
t282 = 0.2e1 * qJ(3,1);
t222 = sin(t282);
t459 = pkin(3) * t222;
t53 = (t262 * (t100 * t253 + t104 * t259) / (t459 / 0.2e1 + t399) + (t193 * t263 - t196 * t264) * (t100 * t259 - t104 * t253) * t413) * t286;
t377 = t53 * t456;
t245 = xDDP(3);
t412 = t134 * t245;
t452 = t239 * pkin(1);
t305 = -t265 / 0.2e1 + t452 / 0.2e1;
t472 = -0.2e1 * t262;
t70 = t108 * t264 + t111 * t263 + t189 * t472;
t44 = t134 * t305 * t70 - t377;
t442 = t70 / 0.2e1;
t42 = -t189 * t412 + (-t44 + t377) * t135 * t442 + (t108 * t247 + t111 * t246) * t134 / 0.2e1;
t467 = t42 * pkin(1);
t50 = t53 ^ 2;
t505 = t50 * pkin(5) + t239 * (pkin(1) * t50 - t96) + t240 * (-0.2e1 * t467 - t95) - t42 * t283;
t212 = qJ(1,2) + pkin(7);
t162 = t243 + t212;
t163 = -t243 + t212;
t107 = -sin(t162) - sin(t163);
t110 = cos(t163) - cos(t162);
t256 = cos(qJ(3,2));
t198 = t256 * pkin(3);
t176 = t198 + pkin(2);
t141 = t190 + t176;
t131 = 0.1e1 / t141;
t132 = 0.1e1 / t141 ^ 2;
t188 = cos(t212);
t250 = sin(qJ(3,2));
t457 = pkin(3) * t250;
t407 = t239 * t176;
t103 = t168 + t407;
t402 = t250 * t166;
t233 = 0.1e1 / t250;
t416 = t131 * t233;
t279 = 0.2e1 * qJ(3,2);
t219 = sin(t279);
t460 = pkin(3) * t219;
t99 = t240 * t176 + t156;
t52 = (t262 * (t103 * t257 + t99 * t251) / (t460 / 0.2e1 + t402) + (t192 * t263 - t195 * t264) * (-t103 * t251 + t99 * t257) * t416) * t286;
t378 = t52 * t457;
t415 = t131 * t245;
t71 = t107 * t264 + t110 * t263 + t188 * t472;
t441 = t71 / 0.2e1;
t45 = t131 * t305 * t71 - t378;
t41 = -t188 * t415 + (-t45 + t378) * t132 * t441 + (t107 * t247 + t110 * t246) * t131 / 0.2e1;
t468 = t41 * pkin(1);
t49 = t52 ^ 2;
t504 = t49 * pkin(5) + t239 * (pkin(1) * t49 - t94) + t240 * (-0.2e1 * t468 - t93) - t41 * t283;
t211 = qJ(1,3) + pkin(7);
t160 = t242 + t211;
t161 = -t242 + t211;
t106 = -sin(t160) - sin(t161);
t109 = cos(t161) - cos(t160);
t254 = cos(qJ(3,3));
t197 = t254 * pkin(3);
t174 = t197 + pkin(2);
t140 = t190 + t174;
t128 = 0.1e1 / t140;
t129 = 0.1e1 / t140 ^ 2;
t187 = cos(t211);
t248 = sin(qJ(3,3));
t458 = pkin(3) * t248;
t408 = t239 * t174;
t102 = t168 + t408;
t232 = 0.1e1 / t248;
t419 = t128 * t232;
t339 = t286 * t419;
t405 = t248 * t166;
t276 = 0.2e1 * qJ(3,3);
t216 = sin(t276);
t461 = pkin(3) * t216;
t98 = t240 * t174 + t156;
t51 = t262 * (t102 * t255 + t98 * t249) * t286 / (t461 / 0.2e1 + t405) + (t191 * t263 - t194 * t264) * (-t102 * t249 + t98 * t255) * t339;
t379 = t51 * t458;
t418 = t128 * t245;
t69 = t106 * t264 + t109 * t263 + t187 * t472;
t43 = t128 * t305 * t69 - t379;
t443 = t69 / 0.2e1;
t40 = -t187 * t418 + (-t43 + t379) * t129 * t443 + (t106 * t247 + t109 * t246) * t128 / 0.2e1;
t469 = t40 * pkin(1);
t48 = t51 ^ 2;
t503 = t48 * pkin(5) + t239 * (pkin(1) * t48 - t92) + t240 * (-0.2e1 * t469 - t91) - t40 * t283;
t502 = 0.4e1 * pkin(2);
t387 = 0.2e1 * pkin(3);
t471 = 0.2e1 * t265;
t424 = t111 * t134;
t501 = t424 / 0.2e1;
t426 = t110 * t131;
t500 = t426 / 0.2e1;
t428 = t109 * t128;
t499 = t428 / 0.2e1;
t430 = t108 * t134;
t498 = t430 / 0.2e1;
t432 = t107 * t131;
t497 = t432 / 0.2e1;
t434 = t106 * t128;
t496 = t434 / 0.2e1;
t73 = (t174 * t255 - t403) * t240 + t125 - t174 * t404;
t495 = t128 * t73;
t75 = (t176 * t257 - t400) * t240 + t126 - t176 * t401;
t494 = t131 * t75;
t77 = (t178 * t259 - t397) * t240 + t127 - t178 * t398;
t493 = t134 * t77;
t492 = t129 * t254;
t491 = t132 * t256;
t490 = t135 * t258;
t287 = pkin(2) ^ 2;
t288 = pkin(1) ^ 2;
t390 = 0.2e1 * t287 + t288;
t273 = 0.2e1 * pkin(7);
t204 = cos(t273);
t182 = t288 * t204;
t285 = pkin(3) ^ 2;
t483 = -t182 - t285 / 0.2e1;
t482 = 0.2e1 * pkin(1);
t480 = -0.4e1 * pkin(3);
t380 = pkin(2) * t190;
t118 = 0.4e1 * t380 + t182 + t390;
t159 = -t265 + t452;
t205 = pkin(7) + qJ(3,3);
t170 = 0.2e1 * t205;
t206 = -pkin(7) + qJ(3,3);
t171 = 0.2e1 * t206;
t173 = t197 + t283;
t183 = sin(t205);
t184 = sin(t206);
t275 = 0.3e1 * qJ(3,3);
t207 = t275 + pkin(7);
t208 = t275 - pkin(7);
t209 = t276 + pkin(7);
t210 = t276 - pkin(7);
t274 = 0.4e1 * qJ(3,3);
t214 = sin(t274);
t215 = sin(t275);
t224 = cos(t275);
t225 = cos(t276);
t271 = 0.2e1 * t288;
t284 = pkin(3) * t285;
t203 = sin(t273);
t238 = t265 ^ 2;
t267 = 0.3e1 * t287;
t371 = pkin(1) * t169;
t391 = t238 + t288;
t296 = -0.8e1 * (0.3e1 / 0.4e1 * t285 + t267 + t391) * t190 - 0.8e1 * (-t371 + 0.3e1 / 0.8e1 * t285 + t288 + t287 / 0.2e1 + t238 / 0.2e1) * t283 + 0.8e1 * (-pkin(2) * t204 + t203 * t265) * t288;
t315 = t483 - t390;
t316 = -0.2e1 * t182 - 0.4e1 * t287 - 0.2e1 * t288 - t285;
t370 = pkin(1) * t168;
t392 = t288 * t203;
t330 = -0.2e1 * t370 + t392;
t332 = pkin(3) * (t224 - t254);
t336 = cos(t209) + cos(t210) - 0.2e1 * t240;
t374 = 0.2e1 * t392;
t375 = -0.6e1 * t166 * t285;
t376 = -0.2e1 * t159 * t285;
t157 = -0.2e1 * t371;
t389 = (0.6e1 * t380 + t157 + t271 + t267 + t238 - t483) * t480;
t395 = t285 * cos(t274);
t451 = (t187 * t471 + sin(t211) * t283 + t249 * t482 + (sin(qJ(1,3) - t206) + sin(qJ(1,3) + t205)) * pkin(3)) / (t248 * t283 + t461 + (t183 + t184) * pkin(1));
t297 = 0.2e1 * t51;
t63 = t128 * t443;
t46 = -t297 + t63;
t47 = t297 + t63;
t90 = t159 * t283 + t330;
t12 = -(t43 * t502 + ((-t51 + t63) * sin(t171) - (t51 + t63) * sin(t170)) * t288 + (-0.2e1 * (t285 + t390) * t216 + pkin(2) * t215 * t480 - t285 * t214) * t51 + (t374 + (t225 * t502 + 0.2e1 * t332) * t265) * t63 + (((-sin(t207) - t184) * t47 + (sin(t208) + t183) * t46) * pkin(3) + (sin(t210) * t46 - sin(t209) * t47) * t283 + t336 * t63 * t471) * pkin(1)) / ((t271 + 0.4e1 * t287) * t225 + t395 + (cos(t171) + cos(t170)) * t288 + t332 * t502 + (t336 * t502 + (cos(t208) + cos(t207) - cos(t206) - cos(t205)) * t387) * pkin(1) + t316) * t51 + (t191 * t246 - t194 * t247) * t73 * t339 + (t245 * t451 - ((t224 * t376 + (t173 * t159 - t90 * t225 + t330) * t387) * t51 + (-t214 * t284 + t215 * t375 + t216 * t389 + t248 * t296) * t63) / (t118 * t225 + t395 / 0.2e1 - 0.2e1 * t173 * t190 + (-pkin(2) * t254 + t166 * t224) * t387 + t315) * t63 / 0.2e1) * t286;
t329 = t51 * t63;
t463 = t239 / 0.2e1;
t470 = pkin(5) / 0.2e1;
t479 = -0.2e1 * pkin(2) * t329 - 0.2e1 * t12 * t470 - 0.2e1 * (t12 * t463 + t240 * t329) * pkin(1);
t177 = t199 + t283;
t280 = 0.4e1 * qJ(3,1);
t220 = sin(t280);
t281 = 0.3e1 * qJ(3,1);
t221 = sin(t281);
t230 = cos(t281);
t298 = pkin(3) * (-t258 * pkin(2) + t166 * t230);
t303 = -0.4e1 * t370 + t374;
t331 = t159 * t387;
t351 = t134 * t442;
t352 = t77 * t413;
t372 = t177 * t190;
t386 = 0.4e1 * pkin(3);
t393 = t285 * cos(t280);
t231 = cos(t282);
t421 = t118 * t231;
t444 = t231 * t90;
t382 = pkin(7) + qJ(3,1);
t384 = -pkin(7) + qJ(3,1);
t449 = (t189 * t471 + sin(t213) * t283 + t253 * t482 + (sin(qJ(1,1) - t384) + sin(qJ(1,1) + t382)) * pkin(3)) / (t252 * t283 + t459 + (sin(t382) + sin(t384)) * pkin(1));
t476 = -0.2e1 * t285 - 0.2e1 * t118;
t18 = -(t44 * t502 + (t222 * t476 - t285 * t220 + (-t166 * t221 - t252 * t190) * t386) * t53 + (-0.2e1 * t444 + (-t230 + t258) * t331 + t303) * t351) / (0.4e1 * t298 + t316 - 0.4e1 * t372 + t393 + 0.2e1 * t421) * t53 + (t245 * t449 - ((t230 * t376 + (t177 * t159 + t330 - t444) * t387) * t53 + (-t220 * t284 + t221 * t375 + t222 * t389 + t252 * t296) * t351) / (t421 + t393 / 0.2e1 - 0.2e1 * t372 + 0.2e1 * t298 + t315) * t351 / 0.2e1 + (t193 * t246 - t196 * t247) * t352) * t286;
t325 = t53 * t351;
t478 = -0.2e1 * pkin(2) * t325 - 0.2e1 * t18 * t470 - 0.2e1 * (t18 * t463 + t240 * t325) * pkin(1);
t175 = t198 + t283;
t277 = 0.4e1 * qJ(3,2);
t217 = sin(t277);
t278 = 0.3e1 * qJ(3,2);
t218 = sin(t278);
t227 = cos(t278);
t299 = pkin(3) * (-t256 * pkin(2) + t166 * t227);
t355 = t131 * t441;
t356 = t75 * t416;
t373 = t175 * t190;
t394 = t285 * cos(t277);
t228 = cos(t279);
t422 = t118 * t228;
t445 = t228 * t90;
t381 = pkin(7) + qJ(3,2);
t383 = -pkin(7) + qJ(3,2);
t450 = (t188 * t471 + sin(t212) * t283 + t251 * t482 + (sin(qJ(1,2) - t383) + sin(qJ(1,2) + t381)) * pkin(3)) / (t250 * t283 + t460 + (sin(t381) + sin(t383)) * pkin(1));
t17 = -(t45 * t502 + (t219 * t476 - t285 * t217 + (-t166 * t218 - t250 * t190) * t386) * t52 + (-0.2e1 * t445 + (-t227 + t256) * t331 + t303) * t355) / (0.4e1 * t299 + t316 - 0.4e1 * t373 + t394 + 0.2e1 * t422) * t52 + (t245 * t450 - ((t227 * t376 + (t175 * t159 + t330 - t445) * t387) * t52 + (-t217 * t284 + t218 * t375 + t219 * t389 + t250 * t296) * t355) / (t422 + t394 / 0.2e1 - 0.2e1 * t373 + 0.2e1 * t299 + t315) * t355 / 0.2e1 + (t192 * t246 - t195 * t247) * t356) * t286;
t327 = t52 * t355;
t477 = -0.2e1 * pkin(2) * t327 - 0.2e1 * t17 * t470 - 0.2e1 * (t17 * t463 + t240 * t327) * pkin(1);
t475 = 0.2e1 * t254;
t474 = 0.2e1 * t256;
t473 = 0.2e1 * t258;
t466 = t51 * pkin(3);
t465 = t52 * pkin(3);
t464 = t53 * pkin(3);
t462 = pkin(1) / 0.2e1;
t66 = t69 ^ 2;
t448 = t129 * t66;
t68 = t71 ^ 2;
t447 = t132 * t68;
t67 = t70 ^ 2;
t446 = t135 * t67;
t440 = t286 * t73;
t439 = t286 * t75;
t438 = t286 * t77;
t437 = t40 * t248;
t436 = t41 * t250;
t435 = t42 * t252;
t433 = t106 / 0.2e1;
t431 = t107 / 0.2e1;
t429 = t108 / 0.2e1;
t427 = t109 / 0.2e1;
t425 = t110 / 0.2e1;
t423 = t111 / 0.2e1;
t420 = t128 * t187;
t417 = t131 * t188;
t414 = t134 * t189;
t411 = t232 * t254;
t410 = t233 * t256;
t409 = t234 * t258;
t396 = t286 / 0.4e1;
t385 = t166 * t387;
t235 = t254 ^ 2;
t369 = pkin(3) * (-t255 * t240 + t404) * t235;
t236 = t256 ^ 2;
t368 = pkin(3) * (-t257 * t240 + t401) * t236;
t237 = t258 ^ 2;
t367 = pkin(3) * (-t259 * t240 + t398) * t237;
t366 = t40 * t495;
t365 = t41 * t494;
t364 = t42 * t493;
t363 = t254 * t451;
t362 = t256 * t450;
t361 = t258 * t449;
t57 = t191 * t369 + (-t508 * t191 + t194 * t458) * t254 + t194 * t405;
t60 = -t194 * t369 + (t191 * t458 + t508 * t194) * t254 + t191 * t405;
t72 = (t174 * t249 + t255 * t265) * t240 + t156 * t249 + t255 * t408;
t97 = t157 + t287 + 0.2e1 * t380 + t391;
t28 = -t191 * g(1) - t194 * g(2) + (-t254 * t72 * t418 + (-t159 * t379 + (t235 * t285 + t254 * t385 + t97) * t63) * t443 * t492) * t232 + (t60 * t247 + t57 * t246 - ((t159 * t248 * t63 - t466) * t254 - t51 * t166) * t466) * t419;
t360 = t28 * t419;
t359 = t73 * t419;
t358 = t448 / 0.4e1;
t58 = t192 * t368 + (-t507 * t192 + t195 * t457) * t256 + t195 * t402;
t61 = -t195 * t368 + (t192 * t457 + t507 * t195) * t256 + t192 * t402;
t74 = (t176 * t251 + t257 * t265) * t240 + t156 * t251 + t257 * t407;
t29 = -t192 * g(1) - t195 * g(2) + (-t256 * t74 * t415 + (-t159 * t378 + (t236 * t285 + t256 * t385 + t97) * t355) * t441 * t491) * t233 + (t61 * t247 + t58 * t246 - ((t159 * t250 * t355 - t465) * t256 - t52 * t166) * t465) * t416;
t357 = t29 * t416;
t354 = t447 / 0.4e1;
t59 = t193 * t367 + (-t506 * t193 + t196 * t456) * t258 + t196 * t399;
t62 = -t196 * t367 + (t193 * t456 + t506 * t196) * t258 + t193 * t399;
t76 = (t178 * t253 + t259 * t265) * t240 + t156 * t253 + t259 * t406;
t30 = -t193 * g(1) - t196 * g(2) + (-t258 * t76 * t412 + (-t159 * t377 + (t237 * t285 + t258 * t385 + t97) * t351) * t442 * t490) * t234 + (t62 * t247 + t59 * t246 - ((t159 * t252 * t351 - t464) * t258 - t53 * t166) * t464) * t413;
t353 = t30 * t413;
t350 = t446 / 0.4e1;
t349 = t191 * t440;
t348 = t192 * t439;
t347 = t193 * t438;
t346 = t194 * t440;
t345 = t195 * t439;
t344 = t196 * t438;
t343 = t72 * t411;
t342 = t74 * t410;
t341 = t76 * t409;
t340 = t128 * t411;
t338 = t131 * t410;
t337 = t134 * t409;
t328 = t66 * t492 * t495;
t326 = t68 * t491 * t494;
t324 = t67 * t490 * t493;
t323 = t191 * t359;
t322 = t194 * t359;
t320 = t192 * t356;
t319 = t195 * t356;
t318 = t193 * t352;
t317 = t196 * t352;
t308 = t40 * t73 * t340;
t307 = t41 * t75 * t338;
t306 = t42 * t77 * t337;
t56 = (-0.2e1 * t237 + 0.1e1) * t350;
t55 = (-0.2e1 * t236 + 0.1e1) * t354;
t54 = (-0.2e1 * t235 + 0.1e1) * t358;
t39 = t467 + t95;
t38 = t468 + t93;
t37 = t469 + t91;
t36 = (t325 * t473 + t435) * t252;
t35 = (t327 * t474 + t436) * t250;
t34 = (t329 * t475 + t437) * t248;
t33 = t435 * t473 + (0.4e1 * t237 - 0.2e1) * t325;
t32 = t436 * t474 + (0.4e1 * t236 - 0.2e1) * t327;
t31 = t437 * t475 + (0.4e1 * t235 - 0.2e1) * t329;
t27 = (pkin(1) * t354 + t94) * t240 - t38 * t239 + pkin(2) * t354 - t41 * pkin(5);
t26 = (pkin(1) * t350 + t96) * t240 - t39 * t239 + pkin(2) * t350 - t42 * pkin(5);
t25 = (pkin(1) * t358 + t92) * t240 - t37 * t239 + pkin(2) * t358 - t40 * pkin(5);
t24 = -t250 * t29 + t27 * t256;
t23 = t27 * t250 + t29 * t256;
t22 = -t252 * t30 + t26 * t258;
t21 = t26 * t252 + t30 * t258;
t20 = -t248 * t28 + t25 * t254;
t19 = t25 * t248 + t28 * t254;
t16 = t18 * t258 - t50 * t252;
t15 = t18 * t252 + t50 * t258;
t14 = t17 * t256 - t49 * t250;
t13 = t17 * t250 + t49 * t256;
t11 = t12 * t254 - t48 * t248;
t10 = t12 * t248 + t48 * t254;
t7 = t252 * t505 + t258 * t478;
t6 = t252 * t478 - t258 * t505;
t5 = t250 * t504 + t256 * t477;
t4 = t250 * t477 - t256 * t504;
t2 = t248 * t503 + t254 * t479;
t1 = t248 * t479 - t254 * t503;
t3 = [t40 * t496 + t41 * t497 + t42 * t498, t91 * t496 + t93 * t497 + t95 * t498, t92 * t496 + t94 * t497 + t96 * t498, t60 * t360 + t61 * t357 + t62 * t353 + (t37 * t434 + t38 * t432 + t39 * t430) * t462, (t194 * t328 + t195 * t326 + t196 * t324) * t396 + t34 * t496 + t35 * t497 + t36 * t498, (-t317 * t56 - t319 * t55 - t322 * t54) * t286 + t31 * t496 + t32 * t497 + t33 * t498, (-t194 * t366 - t195 * t365 - t196 * t364) * t286 + t10 * t496 + t13 * t497 + t15 * t498, (-t194 * t308 - t195 * t307 - t196 * t306) * t286 + t11 * t496 + t14 * t497 + t16 * t498, (-t12 * t322 - t17 * t319 - t18 * t317) * t286, (t6 * t429 + (t16 * t62 - t21 * t344) * t234) * t134 + (t4 * t431 + (t14 * t61 - t23 * t345) * t233) * t131 + (t1 * t433 + (t11 * t60 - t19 * t346) * t232) * t128, (t7 * t429 + (-t15 * t62 - t22 * t344) * t234) * t134 + (t5 * t431 + (-t13 * t61 - t24 * t345) * t233) * t131 + (t2 * t433 + (-t10 * t60 - t20 * t346) * t232) * t128, t247 - g(1); t40 * t499 + t41 * t500 + t42 * t501, t91 * t499 + t93 * t500 + t95 * t501, t92 * t499 + t94 * t500 + t96 * t501, t57 * t360 + t58 * t357 + t59 * t353 + (t37 * t428 + t38 * t426 + t39 * t424) * t462, (-t191 * t328 - t192 * t326 - t193 * t324) * t396 + t34 * t499 + t35 * t500 + t36 * t501, (t318 * t56 + t320 * t55 + t323 * t54) * t286 + t31 * t499 + t32 * t500 + t33 * t501, (t191 * t366 + t192 * t365 + t193 * t364) * t286 + t10 * t499 + t13 * t500 + t15 * t501, (t191 * t308 + t192 * t307 + t193 * t306) * t286 + t11 * t499 + t14 * t500 + t16 * t501, (t12 * t323 + t17 * t320 + t18 * t318) * t286, (t6 * t423 + (t16 * t59 + t21 * t347) * t234) * t134 + (t4 * t425 + (t14 * t58 + t23 * t348) * t233) * t131 + (t1 * t427 + (t11 * t57 + t19 * t349) * t232) * t128, (t7 * t423 + (-t15 * t59 + t22 * t347) * t234) * t134 + (t5 * t425 + (-t13 * t58 + t24 * t348) * t233) * t131 + (t2 * t427 + (-t10 * t57 + t20 * t349) * t232) * t128, t246 - g(2); -t40 * t420 - t41 * t417 - t42 * t414, -t95 * t414 - t93 * t417 - t91 * t420, -t96 * t414 - t94 * t417 - t92 * t420, -t72 * t28 * t340 - t74 * t29 * t338 - t76 * t30 * t337 + (-t37 * t420 - t38 * t417 - t39 * t414) * pkin(1), -t34 * t420 - t35 * t417 - t36 * t414 + (-t248 * t363 * t448 - t250 * t362 * t447 - t252 * t361 * t446) * t396, -t31 * t420 - t32 * t417 - t33 * t414 + (t56 * t449 + t55 * t450 + t54 * t451) * t286, -t10 * t420 - t13 * t417 - t15 * t414 + (t435 * t449 + t436 * t450 + t437 * t451) * t286, -t11 * t420 - t14 * t417 - t16 * t414 + (t361 * t42 + t362 * t41 + t363 * t40) * t286, (t12 * t451 + t17 * t450 + t18 * t449) * t286, (-t16 * t341 - t189 * t6) * t134 + (-t14 * t342 - t188 * t4) * t131 + (-t187 * t1 - t11 * t343) * t128 + (t19 * t451 + t21 * t449 + t23 * t450) * t286, (t15 * t341 - t189 * t7) * t134 + (t13 * t342 - t188 * t5) * t131 + (t10 * t343 - t187 * t2) * t128 + (t20 * t451 + t22 * t449 + t24 * t450) * t286, t245 - g(3);];
tauX_reg  = t3;
