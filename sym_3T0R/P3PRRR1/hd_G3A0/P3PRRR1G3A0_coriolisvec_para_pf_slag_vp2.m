% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRR1G3A0
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
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR1G3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR1G3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR1G3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:50
% EndTime: 2020-03-09 21:06:51
% DurationCPUTime: 1.16s
% Computational Cost: add. (13542->131), mult. (6396->276), div. (2988->5), fcn. (7728->24), ass. (0->147)
t463 = 0.1e1 / pkin(2);
t447 = pkin(7) + qJ(2,1);
t438 = qJ(3,1) + t447;
t426 = sin(t438);
t429 = cos(t438);
t432 = sin(t447);
t435 = cos(t447);
t521 = 0.1e1 / (t426 * t435 - t432 * t429);
t524 = t463 * t521;
t448 = legFrame(3,2);
t442 = cos(t448);
t445 = pkin(7) + qJ(2,3);
t436 = qJ(3,3) + t445;
t424 = sin(t436);
t427 = cos(t436);
t430 = sin(t445);
t433 = cos(t445);
t523 = 0.1e1 / (t424 * t433 - t430 * t427);
t527 = t442 * t523;
t449 = legFrame(2,2);
t443 = cos(t449);
t446 = pkin(7) + qJ(2,2);
t437 = qJ(3,2) + t446;
t425 = sin(t437);
t428 = cos(t437);
t431 = sin(t446);
t434 = cos(t446);
t522 = 0.1e1 / (t425 * t434 - t431 * t428);
t526 = t443 * t522;
t450 = legFrame(1,2);
t444 = cos(t450);
t525 = t444 * t521;
t507 = t522 * t463;
t508 = t523 * t463;
t493 = t427 * t527;
t459 = xDP(1);
t500 = t459 * t463;
t400 = t493 * t500;
t457 = xDP(3);
t461 = 0.1e1 / pkin(3);
t415 = pkin(2) * t433 + pkin(3) * t427;
t488 = t415 * t527;
t482 = t459 * t488;
t475 = t461 * t482;
t499 = t461 * t415;
t439 = sin(t448);
t458 = xDP(2);
t503 = t439 * t458;
t412 = pkin(2) * t430 + pkin(3) * t424;
t506 = t412 * t461;
t388 = t400 + (-t475 + ((-t424 + t506) * t457 + (-t427 + t499) * t503) * t523) * t463;
t496 = t461 * t463;
t391 = (-t482 + (t412 * t457 + t415 * t503) * t523) * t496;
t520 = t388 * t391;
t491 = t428 * t526;
t401 = t491 * t500;
t416 = pkin(2) * t434 + pkin(3) * t428;
t487 = t416 * t526;
t481 = t459 * t487;
t474 = t461 * t481;
t498 = t461 * t416;
t440 = sin(t449);
t502 = t440 * t458;
t413 = pkin(2) * t431 + pkin(3) * t425;
t505 = t413 * t461;
t389 = t401 + (-t474 + ((-t425 + t505) * t457 + (-t428 + t498) * t502) * t522) * t463;
t392 = (-t481 + (t413 * t457 + t416 * t502) * t522) * t496;
t519 = t389 * t392;
t489 = t429 * t525;
t402 = t489 * t500;
t417 = pkin(2) * t435 + pkin(3) * t429;
t486 = t417 * t525;
t480 = t459 * t486;
t473 = t461 * t480;
t497 = t461 * t417;
t441 = sin(t450);
t501 = t441 * t458;
t414 = pkin(2) * t432 + pkin(3) * t426;
t504 = t414 * t461;
t390 = t402 + (-t473 + ((-t426 + t504) * t457 + (-t429 + t497) * t501) * t521) * t463;
t393 = (-t480 + (t414 * t457 + t417 * t501) * t521) * t496;
t518 = t390 * t393;
t397 = t400 + (-t424 * t457 - t427 * t503) * t508;
t394 = t397 ^ 2;
t451 = sin(qJ(3,3));
t454 = cos(qJ(3,3));
t421 = t451 * mrSges(3,1) + t454 * mrSges(3,2);
t517 = t394 * t421;
t398 = t401 + (-t425 * t457 - t428 * t502) * t507;
t395 = t398 ^ 2;
t452 = sin(qJ(3,2));
t455 = cos(qJ(3,2));
t422 = t452 * mrSges(3,1) + t455 * mrSges(3,2);
t516 = t395 * t422;
t399 = t402 + (-t426 * t457 - t429 * t501) * t524;
t396 = t399 ^ 2;
t453 = sin(qJ(3,1));
t456 = cos(qJ(3,1));
t423 = t453 * mrSges(3,1) + t456 * mrSges(3,2);
t515 = t396 * t423;
t514 = t523 * t421;
t513 = t523 * t439;
t512 = t522 * t422;
t511 = t522 * t440;
t510 = t521 * t423;
t509 = t521 * t441;
t495 = 0.2e1 * pkin(2) * pkin(3);
t494 = t415 * t513;
t492 = t416 * t511;
t490 = t417 * t509;
t385 = t400 + (-t475 / 0.2e1 + ((-t424 + t506 / 0.2e1) * t457 + (-t427 + t499 / 0.2e1) * t503) * t523) * t463;
t485 = t385 * t391 * t514;
t386 = t401 + (-t474 / 0.2e1 + ((-t425 + t505 / 0.2e1) * t457 + (-t428 + t498 / 0.2e1) * t502) * t522) * t463;
t484 = t386 * t392 * t512;
t387 = t402 + (-t473 / 0.2e1 + ((-t426 + t504 / 0.2e1) * t457 + (-t429 + t497 / 0.2e1) * t501) * t521) * t463;
t483 = t387 * t393 * t510;
t462 = pkin(2) ^ 2;
t479 = -m(3) * t462 - Ifges(2,3) - Ifges(3,3);
t478 = t427 * t485;
t477 = t428 * t484;
t476 = t429 * t483;
t472 = t424 * t430 + t427 * t433;
t471 = t425 * t431 + t428 * t434;
t470 = t426 * t432 + t429 * t435;
t469 = (-mrSges(3,1) * t454 + mrSges(3,2) * t451) * pkin(2);
t468 = (-mrSges(3,1) * t455 + mrSges(3,2) * t452) * pkin(2);
t467 = (-mrSges(3,1) * t456 + mrSges(3,2) * t453) * pkin(2);
t466 = t472 * pkin(2);
t465 = t471 * pkin(2);
t464 = t470 * pkin(2);
t460 = pkin(3) ^ 2;
t420 = -Ifges(3,3) + t467;
t419 = -Ifges(3,3) + t468;
t418 = -Ifges(3,3) + t469;
t384 = (pkin(3) * t518 + (pkin(3) * t390 + t399 * t464) * t399) * t524;
t383 = (pkin(3) * t519 + (t389 * pkin(3) + t398 * t465) * t398) * t507;
t382 = (pkin(3) * t520 + (t388 * pkin(3) + t397 * t466) * t397) * t508;
t381 = (-(-t470 * t387 * t495 - t390 * t460 - t462 * t399) * t461 * t399 + (pkin(3) + t464) * t518) * t524;
t380 = ((-t471 * t386 * t495 - t389 * t460 - t462 * t398) * t461 * t398 - (pkin(3) + t465) * t519) * t507;
t379 = ((-t472 * t385 * t495 - t388 * t460 - t462 * t397) * t461 * t397 - (pkin(3) + t466) * t520) * t508;
t378 = -Ifges(3,3) * t381 - t420 * t384;
t377 = Ifges(3,3) * t380 - t419 * t383;
t376 = Ifges(3,3) * t379 - t418 * t382;
t375 = -(0.2e1 * t467 + t479) * t384 + t420 * t381;
t374 = -(0.2e1 * t468 + t479) * t383 - t419 * t380;
t373 = -(0.2e1 * t469 + t479) * t382 - t418 * t379;
t1 = [-0.2e1 * t442 * t478 - 0.2e1 * t443 * t477 - 0.2e1 * t444 * t476 + (t373 * t493 + t374 * t491 + t375 * t489) * t463 + (-t488 * t517 - t487 * t516 - t486 * t515 + (-t376 * t488 - t377 * t487 - t378 * t486) * t463) * t461; 0.2e1 * t439 * t478 + 0.2e1 * t440 * t477 + 0.2e1 * t441 * t476 + (-t427 * t373 * t513 - t428 * t374 * t511 - t429 * t375 * t509) * t463 + (t494 * t517 + t492 * t516 + t490 * t515 + (t376 * t494 + t377 * t492 + t378 * t490) * t463) * t461; 0.2e1 * t424 * t485 + 0.2e1 * t425 * t484 + 0.2e1 * t426 * t483 + (-t373 * t424 * t523 - t374 * t425 * t522 - t375 * t426 * t521) * t463 + (t394 * t412 * t514 + t395 * t413 * t512 + t396 * t414 * t510 + (t376 * t412 * t523 + t377 * t413 * t522 + t378 * t414 * t521) * t463) * t461;];
taucX  = t1;
