% Calculate Gravitation load for parallel robot
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR2G3P1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:04
% EndTime: 2020-03-09 21:20:04
% DurationCPUTime: 0.20s
% Computational Cost: add. (219->58), mult. (293->88), div. (30->5), fcn. (210->33), ass. (0->46)
t391 = legFrame(3,2);
t378 = sin(t391);
t381 = cos(t391);
t360 = g(1) * t378 + g(2) * t381;
t363 = g(1) * t381 - g(2) * t378;
t388 = qJ(2,3) + qJ(3,3);
t351 = (mrSges(3,1) * t363 + mrSges(3,2) * t360) * cos(t388) + sin(t388) * (mrSges(3,1) * t360 - mrSges(3,2) * t363);
t385 = 0.1e1 / sin(qJ(3,3));
t401 = t351 * t385;
t392 = legFrame(2,2);
t379 = sin(t392);
t382 = cos(t392);
t361 = g(1) * t379 + g(2) * t382;
t364 = g(1) * t382 - g(2) * t379;
t389 = qJ(2,2) + qJ(3,2);
t352 = (mrSges(3,1) * t364 + mrSges(3,2) * t361) * cos(t389) + sin(t389) * (mrSges(3,1) * t361 - mrSges(3,2) * t364);
t386 = 0.1e1 / sin(qJ(3,2));
t400 = t352 * t386;
t393 = legFrame(1,2);
t380 = sin(t393);
t383 = cos(t393);
t362 = g(1) * t380 + g(2) * t383;
t365 = g(1) * t383 - g(2) * t380;
t390 = qJ(2,1) + qJ(3,1);
t353 = (mrSges(3,1) * t365 + mrSges(3,2) * t362) * cos(t390) + sin(t390) * (mrSges(3,1) * t362 - mrSges(3,2) * t365);
t387 = 0.1e1 / sin(qJ(3,1));
t399 = t353 * t387;
t384 = m(3) * pkin(1) + mrSges(2,1);
t398 = t385 * ((mrSges(2,2) * t360 + t363 * t384) * cos(qJ(2,3)) + sin(qJ(2,3)) * (-mrSges(2,2) * t363 + t360 * t384) + t351);
t397 = t386 * ((mrSges(2,2) * t361 + t364 * t384) * cos(qJ(2,2)) + sin(qJ(2,2)) * (-mrSges(2,2) * t364 + t361 * t384) + t352);
t396 = t387 * ((mrSges(2,2) * t362 + t365 * t384) * cos(qJ(2,1)) + sin(qJ(2,1)) * (-mrSges(2,2) * t365 + t362 * t384) + t353);
t375 = -t391 + qJ(2,3);
t376 = -t392 + qJ(2,2);
t377 = -t393 + qJ(2,1);
t395 = 0.1e1 / pkin(1);
t394 = 0.1e1 / pkin(2);
t374 = qJ(3,1) + t377;
t373 = qJ(3,2) + t376;
t372 = qJ(3,3) + t375;
t371 = cos(t374);
t370 = cos(t373);
t369 = cos(t372);
t368 = sin(t374);
t367 = sin(t373);
t366 = sin(t372);
t1 = [-g(1) * m(4) + (-t366 * t398 - t367 * t397 - t368 * t396 + ((pkin(2) * t368 + pkin(1) * sin(t377)) * t399 + (pkin(2) * t367 + pkin(1) * sin(t376)) * t400 + (pkin(2) * t366 + pkin(1) * sin(t375)) * t401) * t394) * t395; -g(2) * m(4) + (t369 * t398 + t370 * t397 + t371 * t396 + ((-pkin(2) * t371 - pkin(1) * cos(t377)) * t399 + (-pkin(2) * t370 - pkin(1) * cos(t376)) * t400 + (-pkin(2) * t369 - pkin(1) * cos(t375)) * t401) * t394) * t395; -(3 * g(3) * (m(1) + m(2) + m(3))) - g(3) * m(4);];
taugX  = t1;
