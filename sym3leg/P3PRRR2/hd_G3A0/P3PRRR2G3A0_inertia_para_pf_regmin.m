% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRR2G3P1A0
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
%   pkin=[a3,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x8]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR2G3P1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3P1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:06
% EndTime: 2020-03-09 21:20:07
% DurationCPUTime: 0.34s
% Computational Cost: add. (704->90), mult. (336->180), div. (408->9), fcn. (486->18), ass. (0->96)
t230 = -legFrame(3,2) + qJ(2,3);
t227 = qJ(3,3) + t230;
t221 = sin(t227);
t215 = pkin(2) * t221 + pkin(1) * sin(t230);
t242 = sin(qJ(3,3));
t233 = 0.1e1 / t242;
t293 = t215 * t233;
t231 = -legFrame(2,2) + qJ(2,2);
t228 = qJ(3,2) + t231;
t222 = sin(t228);
t216 = pkin(2) * t222 + pkin(1) * sin(t231);
t243 = sin(qJ(3,2));
t235 = 0.1e1 / t243;
t292 = t216 * t235;
t232 = -legFrame(1,2) + qJ(2,1);
t229 = qJ(3,1) + t232;
t223 = sin(t229);
t217 = pkin(2) * t223 + pkin(1) * sin(t232);
t244 = sin(qJ(3,1));
t237 = 0.1e1 / t244;
t291 = t217 * t237;
t224 = cos(t227);
t218 = -pkin(2) * t224 - pkin(1) * cos(t230);
t290 = t218 * t233;
t225 = cos(t228);
t219 = -pkin(2) * t225 - pkin(1) * cos(t231);
t289 = t219 * t235;
t226 = cos(t229);
t220 = -pkin(2) * t226 - pkin(1) * cos(t232);
t288 = t220 * t237;
t287 = t221 * t233;
t286 = t222 * t235;
t285 = t223 * t237;
t284 = t224 * t233;
t283 = t225 * t235;
t282 = t226 * t237;
t245 = cos(qJ(3,3));
t281 = t233 * t245;
t249 = 0.1e1 / pkin(1);
t280 = t233 * t249;
t234 = 0.1e1 / t242 ^ 2;
t279 = t234 * t245;
t246 = cos(qJ(3,2));
t278 = t235 * t246;
t277 = t235 * t249;
t236 = 0.1e1 / t243 ^ 2;
t276 = t236 * t246;
t247 = cos(qJ(3,1));
t275 = t237 * t247;
t274 = t237 * t249;
t238 = 0.1e1 / t244 ^ 2;
t273 = t238 * t247;
t248 = 0.1e1 / pkin(2);
t272 = t248 * t249;
t271 = t221 * t281;
t270 = t221 * t279;
t269 = t222 * t278;
t268 = t222 * t276;
t267 = t223 * t275;
t266 = t223 * t273;
t265 = t224 * t281;
t264 = t224 * t279;
t263 = t225 * t278;
t262 = t225 * t276;
t261 = t226 * t275;
t260 = t226 * t273;
t259 = t233 * t272;
t258 = t235 * t272;
t257 = t237 * t272;
t256 = t221 * t280;
t255 = t222 * t277;
t254 = t223 * t274;
t253 = t224 * t280;
t252 = t225 * t277;
t251 = t226 * t274;
t250 = 0.1e1 / pkin(1) ^ 2;
t214 = t220 * t257;
t213 = t219 * t258;
t212 = t218 * t259;
t211 = t217 * t257;
t210 = t216 * t258;
t209 = t215 * t259;
t208 = t214 + t251;
t207 = t213 + t252;
t206 = t212 + t253;
t205 = t211 - t254;
t204 = t210 - t255;
t203 = t209 - t256;
t202 = t214 + 0.2e1 * t251;
t201 = t213 + 0.2e1 * t252;
t200 = t212 + 0.2e1 * t253;
t199 = t211 - 0.2e1 * t254;
t198 = t210 - 0.2e1 * t255;
t197 = t209 - 0.2e1 * t256;
t196 = (-t221 * t224 * t234 - t222 * t225 * t236 - t223 * t226 * t238) * t250;
t1 = [0, (t221 ^ 2 * t234 + t222 ^ 2 * t236 + t223 ^ 2 * t238) * t250, 0, 0, (-t203 * t287 - t204 * t286 - t205 * t285 + (t203 * t293 + t204 * t292 + t205 * t291) * t248) * t249, -t197 * t271 - t198 * t269 - t199 * t267 + (-t215 * t270 - t216 * t268 - t217 * t266) * t272, t221 * t197 + t222 * t198 + t223 * t199 + (t215 * t287 + t216 * t286 + t217 * t285) * t272, 1; 0, t196, 0, 0, (-t206 * t287 - t207 * t286 - t208 * t285 + (t206 * t293 + t207 * t292 + t208 * t291) * t248) * t249, -t200 * t271 - t201 * t269 - t202 * t267 + (t215 * t264 + t216 * t262 + t217 * t260) * t272, t221 * t200 + t222 * t201 + t223 * t202 + (-t215 * t284 - t216 * t283 - t217 * t282) * t272, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, t196, 0, 0, (t203 * t284 + t204 * t283 + t205 * t282 + (t203 * t290 + t204 * t289 + t205 * t288) * t248) * t249, t197 * t265 + t198 * t263 + t199 * t261 + (-t218 * t270 - t219 * t268 - t220 * t266) * t272, -t224 * t197 - t225 * t198 - t226 * t199 + (t218 * t287 + t219 * t286 + t220 * t285) * t272, 0; 0, (t224 ^ 2 * t234 + t225 ^ 2 * t236 + t226 ^ 2 * t238) * t250, 0, 0, (t206 * t284 + t207 * t283 + t208 * t282 + (t206 * t290 + t207 * t289 + t208 * t288) * t248) * t249, t200 * t265 + t201 * t263 + t202 * t261 + (t218 * t264 + t219 * t262 + t220 * t260) * t272, -t224 * t200 - t225 * t201 - t226 * t202 + (-t218 * t284 - t219 * t283 - t220 * t282) * t272, 1; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 3, 0, 0, 0, 0, 0, 0, 1;];
tau_reg  = t1;
