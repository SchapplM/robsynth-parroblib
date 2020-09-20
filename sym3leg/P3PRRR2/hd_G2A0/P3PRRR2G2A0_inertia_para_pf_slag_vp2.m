% Calculate inertia matrix for parallel robot
% P3PRRR2G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRR2G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:35
% EndTime: 2020-03-09 21:21:35
% DurationCPUTime: 0.38s
% Computational Cost: add. (597->97), mult. (1104->200), div. (405->5), fcn. (963->21), ass. (0->111)
t264 = cos(qJ(3,3));
t265 = cos(qJ(2,3));
t258 = sin(qJ(3,3));
t259 = sin(qJ(2,3));
t298 = t259 * t258;
t227 = (pkin(2) * t264 + pkin(1)) * t265 - pkin(2) * t298;
t271 = 0.1e1 / pkin(2);
t319 = t227 * t271;
t266 = cos(qJ(3,2));
t267 = cos(qJ(2,2));
t260 = sin(qJ(3,2));
t261 = sin(qJ(2,2));
t297 = t261 * t260;
t228 = (pkin(2) * t266 + pkin(1)) * t267 - pkin(2) * t297;
t318 = t228 * t271;
t268 = cos(qJ(3,1));
t269 = cos(qJ(2,1));
t262 = sin(qJ(3,1));
t263 = sin(qJ(2,1));
t296 = t263 * t262;
t229 = (pkin(2) * t268 + pkin(1)) * t269 - pkin(2) * t296;
t317 = t229 * t271;
t242 = sin(qJ(2,3) + qJ(3,3));
t239 = t259 * pkin(1) + pkin(2) * t242;
t252 = 0.1e1 / t258;
t316 = t239 * t252;
t315 = t239 * t271;
t243 = sin(qJ(2,2) + qJ(3,2));
t240 = t261 * pkin(1) + pkin(2) * t243;
t253 = 0.1e1 / t260;
t314 = t240 * t253;
t313 = t240 * t271;
t244 = sin(qJ(2,1) + qJ(3,1));
t241 = t263 * pkin(1) + pkin(2) * t244;
t254 = 0.1e1 / t262;
t312 = t241 * t254;
t311 = t241 * t271;
t310 = t242 * t252;
t309 = t243 * t253;
t308 = t244 * t254;
t255 = legFrame(3,2);
t245 = sin(t255);
t307 = t245 * t252;
t256 = legFrame(2,2);
t246 = sin(t256);
t306 = t246 * t253;
t257 = legFrame(1,2);
t247 = sin(t257);
t305 = t247 * t254;
t248 = cos(t255);
t304 = t248 * t252;
t249 = cos(t256);
t303 = t249 * t253;
t250 = cos(t257);
t302 = t250 * t254;
t272 = 0.1e1 / pkin(1);
t301 = t252 * t272;
t300 = t253 * t272;
t299 = t254 * t272;
t295 = m(3) * pkin(1) ^ 2 + Ifges(2,3) + Ifges(3,3);
t294 = t227 * t307;
t293 = t227 * t304;
t292 = t228 * t306;
t291 = t228 * t303;
t290 = t229 * t305;
t289 = t229 * t302;
t233 = t265 * t264 - t298;
t288 = t233 * t307;
t287 = t233 * t304;
t234 = t267 * t266 - t297;
t286 = t234 * t306;
t285 = t234 * t303;
t235 = t269 * t268 - t296;
t284 = t235 * t305;
t283 = t235 * t302;
t282 = (mrSges(3,1) * t264 - mrSges(3,2) * t258) * pkin(1);
t281 = (mrSges(3,1) * t266 - mrSges(3,2) * t260) * pkin(1);
t280 = (mrSges(3,1) * t268 - mrSges(3,2) * t262) * pkin(1);
t236 = Ifges(3,3) + t282;
t279 = (-Ifges(3,3) * t319 + t233 * t236) * t301;
t237 = Ifges(3,3) + t281;
t278 = (-Ifges(3,3) * t318 + t234 * t237) * t300;
t238 = Ifges(3,3) + t280;
t277 = (-Ifges(3,3) * t317 + t235 * t238) * t299;
t230 = 0.2e1 * t282 + t295;
t276 = (t230 * t233 - t236 * t319) * t301;
t231 = 0.2e1 * t281 + t295;
t275 = (t231 * t234 - t237 * t318) * t300;
t232 = 0.2e1 * t280 + t295;
t274 = (t232 * t235 - t238 * t317) * t299;
t251 = m(1) + m(2) + m(3);
t273 = (-t245 * t248 - t246 * t249 - t247 * t250) * t251;
t226 = (Ifges(3,3) * t311 - t238 * t244) * t299;
t225 = (Ifges(3,3) * t313 - t237 * t243) * t300;
t224 = (Ifges(3,3) * t315 - t236 * t242) * t301;
t223 = (-t232 * t244 + t238 * t311) * t299;
t222 = (-t231 * t243 + t237 * t313) * t300;
t221 = (-t230 * t242 + t236 * t315) * t301;
t220 = t250 * t277;
t219 = t249 * t278;
t218 = t248 * t279;
t217 = t247 * t277;
t216 = t246 * t278;
t215 = t245 * t279;
t214 = t250 * t274;
t213 = t249 * t275;
t212 = t248 * t276;
t211 = t247 * t274;
t210 = t246 * t275;
t209 = t245 * t276;
t1 = [m(4) + (t248 ^ 2 + t249 ^ 2 + t250 ^ 2) * t251 + (t209 * t288 + t210 * t286 + t211 * t284 + (-t215 * t294 - t216 * t292 - t217 * t290) * t271) * t272, t273 + (t209 * t287 + t210 * t285 + t211 * t283 + (-t215 * t293 - t216 * t291 - t217 * t289) * t271) * t272, (-t209 * t310 - t210 * t309 - t211 * t308 + (t215 * t316 + t216 * t314 + t217 * t312) * t271) * t272; t273 + (t212 * t288 + t213 * t286 + t214 * t284 + (-t218 * t294 - t219 * t292 - t220 * t290) * t271) * t272, m(4) + (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) * t251 + (t212 * t287 + t213 * t285 + t214 * t283 + (-t218 * t293 - t219 * t291 - t220 * t289) * t271) * t272, (-t212 * t310 - t213 * t309 - t214 * t308 + (t218 * t316 + t219 * t314 + t220 * t312) * t271) * t272; (t221 * t288 + t222 * t286 + t223 * t284 + (-t224 * t294 - t225 * t292 - t226 * t290) * t271) * t272, (t221 * t287 + t222 * t285 + t223 * t283 + (-t224 * t293 - t225 * t291 - t226 * t289) * t271) * t272, m(4) + (-t221 * t310 - t222 * t309 - t223 * t308 + (t224 * t316 + t225 * t314 + t226 * t312) * t271) * t272;];
MX  = t1;
