% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRR1G2P2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2P2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:22
% EndTime: 2020-03-09 21:18:23
% DurationCPUTime: 1.08s
% Computational Cost: add. (13542->113), mult. (6396->233), div. (2988->5), fcn. (7728->24), ass. (0->140)
t481 = 0.1e1 / pkin(2);
t465 = pkin(7) + qJ(2,1);
t456 = qJ(3,1) + t465;
t444 = sin(t456);
t447 = cos(t456);
t450 = sin(t465);
t453 = cos(t465);
t539 = 0.1e1 / (-t453 * t444 + t447 * t450);
t547 = t481 * t539;
t464 = pkin(7) + qJ(2,2);
t455 = qJ(3,2) + t464;
t443 = sin(t455);
t446 = cos(t455);
t449 = sin(t464);
t452 = cos(t464);
t540 = 0.1e1 / (-t452 * t443 + t446 * t449);
t546 = t481 * t540;
t463 = pkin(7) + qJ(2,3);
t454 = qJ(3,3) + t463;
t442 = sin(t454);
t445 = cos(t454);
t448 = sin(t463);
t451 = cos(t463);
t541 = 0.1e1 / (-t451 * t442 + t445 * t448);
t545 = t481 * t541;
t544 = t442 * t541;
t543 = t443 * t540;
t542 = t444 * t539;
t466 = legFrame(3,2);
t460 = cos(t466);
t500 = t460 * t544;
t457 = sin(t466);
t506 = t457 * t544;
t477 = xDP(1);
t518 = t477 * t481;
t476 = xDP(2);
t519 = t476 * t481;
t475 = xDP(3);
t520 = t475 * t481;
t532 = t541 * t445;
t409 = -t500 * t518 + t506 * t519 - t520 * t532;
t430 = pkin(2) * t448 + pkin(3) * t442;
t479 = 0.1e1 / pkin(3);
t517 = t479 * t481;
t505 = t541 * t517;
t523 = t460 * t477;
t526 = t457 * t476;
t433 = pkin(2) * t451 + pkin(3) * t445;
t529 = t433 * t475;
t397 = (t529 / 0.2e1 + (t523 / 0.2e1 - t526 / 0.2e1) * t430) * t505 + t409;
t403 = (t529 + (t523 - t526) * t430) * t505;
t400 = t403 + t409;
t478 = pkin(3) ^ 2;
t480 = pkin(2) ^ 2;
t490 = t442 * t448 + t445 * t451;
t484 = t490 * pkin(2);
t516 = 0.2e1 * pkin(2) * pkin(3);
t535 = t400 * t403;
t391 = ((-t490 * t397 * t516 - t400 * t478 - t480 * t409) * t479 * t409 - (pkin(3) + t484) * t535) * t545;
t396 = (pkin(3) * t535 + (t400 * pkin(3) + t409 * t484) * t409) * t545;
t469 = sin(qJ(3,3));
t472 = cos(qJ(3,3));
t487 = (-mrSges(3,1) * t472 + mrSges(3,2) * t469) * pkin(2);
t436 = -Ifges(3,3) + t487;
t538 = (-Ifges(3,3) * t391 + t436 * t396) * t541;
t467 = legFrame(2,2);
t461 = cos(t467);
t499 = t461 * t543;
t458 = sin(t467);
t504 = t458 * t543;
t531 = t540 * t446;
t410 = -t499 * t518 + t504 * t519 - t520 * t531;
t431 = pkin(2) * t449 + pkin(3) * t443;
t503 = t540 * t517;
t522 = t461 * t477;
t525 = t458 * t476;
t434 = pkin(2) * t452 + pkin(3) * t446;
t528 = t434 * t475;
t398 = (t528 / 0.2e1 + (t522 / 0.2e1 - t525 / 0.2e1) * t431) * t503 + t410;
t404 = (t528 + (t522 - t525) * t431) * t503;
t401 = t404 + t410;
t489 = t443 * t449 + t446 * t452;
t483 = t489 * pkin(2);
t534 = t401 * t404;
t392 = ((-t489 * t398 * t516 - t401 * t478 - t480 * t410) * t479 * t410 - (pkin(3) + t483) * t534) * t546;
t394 = (pkin(3) * t534 + (t401 * pkin(3) + t410 * t483) * t410) * t546;
t470 = sin(qJ(3,2));
t473 = cos(qJ(3,2));
t486 = (-mrSges(3,1) * t473 + mrSges(3,2) * t470) * pkin(2);
t437 = -Ifges(3,3) + t486;
t537 = (-Ifges(3,3) * t392 + t437 * t394) * t540;
t468 = legFrame(1,2);
t462 = cos(t468);
t498 = t462 * t542;
t459 = sin(t468);
t502 = t459 * t542;
t530 = t539 * t447;
t411 = -t498 * t518 + t502 * t519 - t520 * t530;
t432 = pkin(2) * t450 + pkin(3) * t444;
t501 = t539 * t517;
t521 = t462 * t477;
t524 = t459 * t476;
t435 = pkin(2) * t453 + pkin(3) * t447;
t527 = t435 * t475;
t399 = (t527 / 0.2e1 + (t521 / 0.2e1 - t524 / 0.2e1) * t432) * t501 + t411;
t405 = (t527 + (t521 - t524) * t432) * t501;
t402 = t405 + t411;
t488 = t444 * t450 + t447 * t453;
t482 = t488 * pkin(2);
t533 = t402 * t405;
t393 = ((-t488 * t399 * t516 - t402 * t478 - t480 * t411) * t479 * t411 - (pkin(3) + t482) * t533) * t547;
t395 = (pkin(3) * t533 + (pkin(3) * t402 + t411 * t482) * t411) * t547;
t471 = sin(qJ(3,1));
t474 = cos(qJ(3,1));
t485 = (-mrSges(3,1) * t474 + mrSges(3,2) * t471) * pkin(2);
t438 = -Ifges(3,3) + t485;
t536 = (-Ifges(3,3) * t393 + t438 * t395) * t539;
t515 = t430 * t538;
t514 = t431 * t537;
t513 = t432 * t536;
t439 = t469 * mrSges(3,1) + t472 * mrSges(3,2);
t512 = t397 * t403 * t439;
t440 = t470 * mrSges(3,1) + t473 * mrSges(3,2);
t511 = t398 * t404 * t440;
t441 = t471 * mrSges(3,1) + t474 * mrSges(3,2);
t510 = t399 * t405 * t441;
t509 = t409 ^ 2 * t541 * t439;
t508 = t410 ^ 2 * t540 * t440;
t507 = t411 ^ 2 * t539 * t441;
t497 = t430 * t509;
t496 = t431 * t508;
t495 = t432 * t507;
t494 = -m(3) * t480 - Ifges(2,3) - Ifges(3,3);
t493 = 0.2e1 * t541 * t512;
t492 = 0.2e1 * t540 * t511;
t491 = 0.2e1 * t539 * t510;
t387 = (0.2e1 * t485 + t494) * t395 + t438 * t393;
t386 = (0.2e1 * t486 + t494) * t394 + t437 * t392;
t385 = (0.2e1 * t487 + t494) * t396 + t436 * t391;
t1 = [t460 * t442 * t493 + t461 * t443 * t492 + t462 * t444 * t491 + (-t385 * t500 - t386 * t499 - t387 * t498) * t481 + (t460 * t497 + t461 * t496 + t462 * t495 + (t460 * t515 + t461 * t514 + t462 * t513) * t481) * t479; -0.2e1 * t506 * t512 - 0.2e1 * t504 * t511 - 0.2e1 * t502 * t510 + (t385 * t506 + t386 * t504 + t387 * t502) * t481 + (-t457 * t497 - t458 * t496 - t459 * t495 + (-t457 * t515 - t458 * t514 - t459 * t513) * t481) * t479; t445 * t493 + t446 * t492 + t447 * t491 + (-t385 * t532 - t386 * t531 - t387 * t530) * t481 + (t433 * t509 + t434 * t508 + t435 * t507 + (t433 * t538 + t434 * t537 + t435 * t536) * t481) * t479;];
taucX  = t1;
