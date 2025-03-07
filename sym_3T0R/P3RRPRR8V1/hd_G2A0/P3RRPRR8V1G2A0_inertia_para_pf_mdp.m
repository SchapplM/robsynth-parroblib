% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR8V1G2A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-04 17:05
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RRPRR8V1G2A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_mdp: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:05:19
% EndTime: 2022-11-04 17:05:26
% DurationCPUTime: 7.33s
% Computational Cost: add. (6840->501), mult. (12968->1037), div. (1527->18), fcn. (11964->23), ass. (0->410)
t2760 = 2 * MDP(13);
t2562 = cos(pkin(5));
t2536 = pkin(2) * t2562 + pkin(1);
t2566 = legFrame(3,2);
t2540 = sin(t2566);
t2543 = cos(t2566);
t2569 = sin(qJ(2,3));
t2575 = cos(qJ(2,3));
t2570 = sin(qJ(1,3));
t2561 = sin(pkin(5));
t2825 = pkin(2) * t2561;
t2753 = t2570 * t2825;
t2793 = t2536 * t2570;
t2485 = (-t2540 * t2793 + t2543 * t2825) * t2575 + (t2536 * t2543 + t2540 * t2753) * t2569;
t2488 = (t2540 * t2825 + t2543 * t2793) * t2575 + (t2536 * t2540 - t2543 * t2753) * t2569;
t2836 = t2485 * t2488;
t2567 = legFrame(2,2);
t2541 = sin(t2567);
t2544 = cos(t2567);
t2571 = sin(qJ(2,2));
t2577 = cos(qJ(2,2));
t2572 = sin(qJ(1,2));
t2752 = t2572 * t2825;
t2792 = t2536 * t2572;
t2486 = (-t2541 * t2792 + t2544 * t2825) * t2577 + (t2536 * t2544 + t2541 * t2752) * t2571;
t2489 = (t2541 * t2825 + t2544 * t2792) * t2577 + (t2536 * t2541 - t2544 * t2752) * t2571;
t2835 = t2486 * t2489;
t2568 = legFrame(1,2);
t2542 = sin(t2568);
t2545 = cos(t2568);
t2573 = sin(qJ(2,1));
t2579 = cos(qJ(2,1));
t2574 = sin(qJ(1,1));
t2751 = t2574 * t2825;
t2791 = t2536 * t2574;
t2487 = (-t2542 * t2791 + t2545 * t2825) * t2579 + (t2536 * t2545 + t2542 * t2751) * t2573;
t2490 = (t2542 * t2825 + t2545 * t2791) * t2579 + (t2536 * t2542 - t2545 * t2751) * t2573;
t2834 = t2487 * t2490;
t2833 = -0.2e1 * pkin(1);
t2832 = 2 * MDP(5);
t2831 = 2 * MDP(6);
t2830 = 2 * MDP(7);
t2829 = pkin(1) * t2562;
t2828 = pkin(1) * t2575;
t2827 = pkin(1) * t2577;
t2826 = pkin(1) * t2579;
t2824 = 0.2e1 * pkin(1);
t2775 = t2561 * t2569;
t2518 = -pkin(2) * t2775 + t2536 * t2575;
t2509 = 0.1e1 / t2518;
t2823 = t2485 * t2509;
t2563 = pkin(4) + qJ(3,3);
t2546 = 0.1e1 / t2563;
t2822 = t2485 * t2546;
t2774 = t2561 * t2571;
t2519 = -pkin(2) * t2774 + t2536 * t2577;
t2511 = 0.1e1 / t2519;
t2821 = t2486 * t2511;
t2564 = pkin(4) + qJ(3,2);
t2548 = 0.1e1 / t2564;
t2820 = t2486 * t2548;
t2773 = t2561 * t2573;
t2520 = -pkin(2) * t2773 + t2536 * t2579;
t2513 = 0.1e1 / t2520;
t2819 = t2487 * t2513;
t2565 = pkin(4) + qJ(3,1);
t2550 = 0.1e1 / t2565;
t2818 = t2487 * t2550;
t2817 = t2488 * t2509;
t2816 = t2488 * t2546;
t2815 = t2489 * t2511;
t2814 = t2489 * t2548;
t2813 = t2490 * t2513;
t2812 = t2490 * t2550;
t2533 = pkin(2) * cos(qJ(2,3) + pkin(5)) + t2828;
t2576 = cos(qJ(1,3));
t2506 = t2533 * t2576 + t2563 * t2570;
t2811 = (-t2576 * t2828 + t2506) * t2546 ^ 2;
t2534 = pkin(2) * cos(qJ(2,2) + pkin(5)) + t2827;
t2578 = cos(qJ(1,2));
t2507 = t2534 * t2578 + t2564 * t2572;
t2810 = (-t2578 * t2827 + t2507) * t2548 ^ 2;
t2535 = pkin(2) * cos(qJ(2,1) + pkin(5)) + t2826;
t2580 = cos(qJ(1,1));
t2508 = t2535 * t2580 + t2565 * t2574;
t2809 = (-t2580 * t2826 + t2508) * t2550 ^ 2;
t2547 = 0.1e1 / t2563 ^ 2;
t2808 = t2506 * t2547;
t2549 = 0.1e1 / t2564 ^ 2;
t2807 = t2507 * t2549;
t2551 = 0.1e1 / t2565 ^ 2;
t2806 = t2508 * t2551;
t2805 = t2509 * t2546;
t2510 = 0.1e1 / t2518 ^ 2;
t2804 = t2510 * t2547;
t2803 = t2511 * t2548;
t2512 = 0.1e1 / t2519 ^ 2;
t2802 = t2512 * t2549;
t2801 = t2513 * t2550;
t2514 = 0.1e1 / t2520 ^ 2;
t2800 = t2514 * t2551;
t2527 = 0.1e1 / t2533;
t2799 = t2527 * t2540;
t2798 = t2527 * t2543;
t2529 = 0.1e1 / t2534;
t2797 = t2529 * t2541;
t2796 = t2529 * t2544;
t2531 = 0.1e1 / t2535;
t2795 = t2531 * t2542;
t2794 = t2531 * t2545;
t2790 = t2546 * t2561;
t2789 = t2546 * t2562;
t2788 = t2546 * t2576;
t2556 = t2576 ^ 2;
t2787 = t2547 * t2556;
t2786 = t2547 * t2576;
t2785 = t2548 * t2561;
t2784 = t2548 * t2562;
t2783 = t2548 * t2578;
t2558 = t2578 ^ 2;
t2782 = t2549 * t2558;
t2781 = t2549 * t2578;
t2780 = t2550 * t2561;
t2779 = t2550 * t2562;
t2778 = t2550 * t2580;
t2560 = t2580 ^ 2;
t2777 = t2551 * t2560;
t2776 = t2551 * t2580;
t2772 = t2561 * t2575;
t2771 = t2561 * t2577;
t2770 = t2561 * t2579;
t2769 = t2562 * t2569;
t2768 = t2562 * t2571;
t2767 = t2562 * t2573;
t2766 = t2569 * t2575;
t2765 = t2571 * t2577;
t2764 = t2573 * t2579;
t2763 = t2561 * t2833;
t2762 = -0.2e1 * t2829;
t2761 = 0.2e1 * t2829;
t2759 = t2509 * t2828;
t2758 = t2511 * t2827;
t2757 = t2513 * t2826;
t2756 = pkin(1) * t2527 * t2569;
t2755 = pkin(1) * t2529 * t2571;
t2754 = pkin(1) * t2531 * t2573;
t2750 = qJ(3,1) * t2531 * t2561;
t2749 = qJ(3,2) * t2529 * t2561;
t2748 = qJ(3,3) * t2527 * t2561;
t2497 = t2518 * t2570 - t2563 * t2576;
t2515 = pkin(2) * t2772 + t2536 * t2569;
t2473 = t2497 * t2543 + t2515 * t2540;
t2747 = t2473 * t2485 * t2547;
t2474 = -t2497 * t2540 + t2515 * t2543;
t2746 = t2474 * t2488 * t2547;
t2498 = t2519 * t2572 - t2564 * t2578;
t2516 = pkin(2) * t2771 + t2536 * t2571;
t2475 = t2498 * t2544 + t2516 * t2541;
t2745 = t2475 * t2486 * t2549;
t2476 = -t2498 * t2541 + t2516 * t2544;
t2744 = t2476 * t2489 * t2549;
t2499 = t2520 * t2574 - t2565 * t2580;
t2517 = pkin(2) * t2770 + t2536 * t2573;
t2477 = t2499 * t2545 + t2517 * t2542;
t2743 = t2477 * t2487 * t2551;
t2478 = -t2499 * t2542 + t2517 * t2545;
t2742 = t2478 * t2490 * t2551;
t2741 = t2485 * t2805;
t2740 = t2486 * t2803;
t2739 = t2487 * t2801;
t2738 = t2488 * t2805;
t2737 = t2489 * t2803;
t2736 = t2490 * t2801;
t2555 = t2575 ^ 2;
t2584 = pkin(1) ^ 2;
t2690 = qJ(3,3) ^ 2 + t2555 * t2584;
t2491 = (-t2506 * t2828 + t2690 * t2576) * t2546;
t2735 = t2491 * t2805;
t2557 = t2577 ^ 2;
t2689 = qJ(3,2) ^ 2 + t2557 * t2584;
t2492 = (-t2507 * t2827 + t2689 * t2578) * t2548;
t2734 = t2492 * t2803;
t2559 = t2579 ^ 2;
t2688 = qJ(3,1) ^ 2 + t2559 * t2584;
t2493 = (-t2508 * t2826 + t2688 * t2580) * t2550;
t2733 = t2493 * t2801;
t2521 = -t2562 * t2575 + t2775;
t2732 = t2521 * t2808;
t2524 = t2769 + t2772;
t2731 = t2524 * t2808;
t2522 = -t2562 * t2577 + t2774;
t2730 = t2522 * t2807;
t2525 = t2768 + t2771;
t2729 = t2525 * t2807;
t2523 = -t2562 * t2579 + t2773;
t2728 = t2523 * t2806;
t2526 = t2767 + t2770;
t2727 = t2526 * t2806;
t2726 = t2527 * t2805;
t2725 = t2569 * t2805;
t2724 = t2509 * t2786;
t2552 = t2569 ^ 2;
t2723 = t2552 * t2804;
t2722 = t2529 * t2803;
t2721 = t2571 * t2803;
t2720 = t2511 * t2781;
t2553 = t2571 ^ 2;
t2719 = t2553 * t2802;
t2718 = t2531 * t2801;
t2717 = t2573 * t2801;
t2716 = t2513 * t2776;
t2554 = t2573 ^ 2;
t2715 = t2554 * t2800;
t2714 = t2521 * t2786;
t2713 = t2522 * t2781;
t2712 = t2523 * t2776;
t2711 = t2524 * t2786;
t2710 = t2525 * t2781;
t2709 = t2526 * t2776;
t2708 = t2527 * t2788;
t2707 = t2529 * t2783;
t2706 = t2531 * t2778;
t2705 = t2547 * t2766;
t2704 = t2549 * t2765;
t2703 = t2551 * t2764;
t2702 = t2527 * t2763;
t2701 = t2527 * t2761;
t2700 = t2529 * t2763;
t2699 = t2529 * t2761;
t2698 = t2531 * t2763;
t2697 = t2531 * t2761;
t2696 = t2555 * t2763;
t2695 = t2557 * t2763;
t2694 = t2559 * t2763;
t2693 = 0.2e1 * qJ(3,1) * t2801;
t2692 = 0.2e1 * qJ(3,2) * t2803;
t2691 = 0.2e1 * qJ(3,3) * t2805;
t2687 = t2540 * t2756;
t2686 = t2543 * t2756;
t2685 = t2541 * t2755;
t2684 = t2544 * t2755;
t2683 = t2542 * t2754;
t2682 = t2545 * t2754;
t2681 = qJ(3,1) * t2717;
t2680 = qJ(3,2) * t2721;
t2679 = qJ(3,3) * t2725;
t2678 = t2804 * t2836;
t2677 = t2485 * t2724;
t2676 = t2802 * t2835;
t2675 = t2486 * t2720;
t2674 = t2800 * t2834;
t2673 = t2487 * t2716;
t2672 = t2488 * t2724;
t2671 = t2489 * t2720;
t2670 = t2490 * t2716;
t2669 = t2509 * t2732;
t2668 = t2509 * t2731;
t2667 = t2511 * t2730;
t2666 = t2511 * t2729;
t2665 = t2513 * t2728;
t2664 = t2513 * t2727;
t2663 = t2527 * t2725;
t2662 = t2575 * t2726;
t2661 = t2552 * t2724;
t2660 = t2510 * t2705;
t2659 = t2529 * t2721;
t2658 = t2577 * t2722;
t2657 = t2553 * t2720;
t2656 = t2512 * t2704;
t2655 = t2531 * t2717;
t2654 = t2579 * t2718;
t2653 = t2554 * t2716;
t2652 = t2514 * t2703;
t2651 = t2569 * t2708;
t2650 = t2575 * t2708;
t2649 = t2571 * t2707;
t2648 = t2577 * t2707;
t2647 = t2573 * t2706;
t2646 = t2579 * t2706;
t2645 = t2725 * t2833;
t2644 = t2721 * t2833;
t2643 = t2717 * t2833;
t2642 = t2509 * t2690;
t2641 = t2511 * t2689;
t2640 = t2688 * t2513;
t2639 = t2540 * t2663;
t2638 = t2543 * t2663;
t2637 = t2509 * t2576 * t2705;
t2636 = t2541 * t2659;
t2635 = t2544 * t2659;
t2634 = t2511 * t2578 * t2704;
t2633 = t2542 * t2655;
t2632 = t2545 * t2655;
t2631 = t2513 * t2580 * t2703;
t2630 = t2540 * t2651;
t2629 = t2543 * t2651;
t2628 = t2541 * t2649;
t2627 = t2544 * t2649;
t2626 = t2542 * t2647;
t2625 = t2545 * t2647;
t2587 = t2625 + t2627 + t2629;
t2624 = (t2485 * t2637 + t2486 * t2634 + t2487 * t2631) * t2832 + (t2485 * t2661 + t2486 * t2657 + t2487 * t2653) * MDP(4) + (t2673 + t2675 + t2677) * MDP(1) + t2587 * MDP(6) + (t2543 * t2650 + t2544 * t2648 + t2545 * t2646) * MDP(7);
t2588 = t2626 + t2628 + t2630;
t2623 = (t2488 * t2637 + t2489 * t2634 + t2490 * t2631) * t2832 + (t2488 * t2661 + t2489 * t2657 + t2490 * t2653) * MDP(4) + (t2670 + t2671 + t2672) * MDP(1) + t2588 * MDP(6) + (t2540 * t2650 + t2541 * t2648 + t2542 * t2646) * MDP(7);
t2622 = t2696 * t2805;
t2621 = t2695 * t2803;
t2620 = t2694 * t2801;
t2528 = 0.1e1 / t2533 ^ 2;
t2530 = 0.1e1 / t2534 ^ 2;
t2532 = 0.1e1 / t2535 ^ 2;
t2592 = (t2487 * t2542 + t2490 * t2545) * t2718;
t2593 = (t2486 * t2541 + t2489 * t2544) * t2722;
t2594 = (t2485 * t2540 + t2488 * t2543) * t2726;
t2619 = (t2569 * t2594 + t2571 * t2593 + t2573 * t2592) * MDP(6) + (t2575 * t2594 + t2577 * t2593 + t2579 * t2592) * MDP(7) + (t2652 * t2834 + t2656 * t2835 + t2660 * t2836) * t2832 + (t2552 * t2678 + t2553 * t2676 + t2554 * t2674) * MDP(4) + (t2674 + t2676 + t2678) * MDP(1) + (t2528 * t2540 * t2543 + t2530 * t2541 * t2544 + t2532 * t2542 * t2545) * MDP(8);
t2618 = qJ(3,1) * t2526;
t2617 = qJ(3,2) * t2525;
t2616 = qJ(3,3) * t2524;
t2615 = t2547 * (MDP(11) * t2521 + MDP(12) * t2524);
t2614 = t2549 * (MDP(11) * t2522 + MDP(12) * t2525);
t2613 = t2551 * (MDP(11) * t2523 + MDP(12) * t2526);
t2612 = t2531 * t2618;
t2611 = t2550 * t2618;
t2610 = qJ(3,1) * t2550 * t2523;
t2609 = t2529 * t2617;
t2608 = t2548 * t2617;
t2607 = qJ(3,2) * t2548 * t2522;
t2606 = t2527 * t2616;
t2605 = t2546 * t2616;
t2604 = qJ(3,3) * t2546 * t2521;
t2603 = t2513 * t2611;
t2602 = t2513 * t2610;
t2601 = t2511 * t2608;
t2600 = t2511 * t2607;
t2599 = t2509 * t2605;
t2598 = t2509 * t2604;
t2597 = (t2555 * t2562 - t2561 * t2766) * t2824;
t2596 = (t2557 * t2562 - t2561 * t2765) * t2824;
t2595 = (t2559 * t2562 - t2561 * t2764) * t2824;
t2591 = t2509 * t2597;
t2590 = t2511 * t2596;
t2589 = t2513 * t2595;
t2586 = t2485 * t2638 + t2486 * t2635 + t2487 * t2632;
t2585 = t2488 * t2639 + t2489 * t2636 + t2490 * t2633;
t2505 = t2580 * t2611;
t2504 = t2580 * t2610;
t2503 = t2578 * t2608;
t2502 = t2578 * t2607;
t2501 = t2576 * t2605;
t2500 = t2576 * t2604;
t2484 = t2490 ^ 2;
t2483 = t2489 ^ 2;
t2482 = t2488 ^ 2;
t2481 = t2487 ^ 2;
t2480 = t2486 ^ 2;
t2479 = t2485 ^ 2;
t2467 = (t2580 * t2694 + (t2573 * t2580 * t2762 + t2508 * t2561) * t2579 + t2508 * t2767) * t2550;
t2466 = (t2578 * t2695 + (t2571 * t2578 * t2762 + t2507 * t2561) * t2577 + t2507 * t2768) * t2548;
t2465 = (t2576 * t2696 + (t2569 * t2576 * t2762 + t2506 * t2561) * t2575 + t2506 * t2769) * t2546;
t2464 = (t2523 * t2508 + t2580 * t2595) * t2550;
t2463 = (t2522 * t2507 + t2578 * t2596) * t2548;
t2462 = (t2521 * t2506 + t2576 * t2597) * t2546;
t2461 = t2487 * t2693 - t2682;
t2460 = t2486 * t2692 - t2684;
t2459 = t2485 * t2691 - t2686;
t2458 = t2490 * t2693 - t2683;
t2457 = t2489 * t2692 - t2685;
t2456 = t2488 * t2691 - t2687;
t2455 = pkin(1) * t2795 - t2490 * t2681;
t2454 = pkin(1) * t2794 - t2487 * t2681;
t2453 = pkin(1) * t2797 - t2489 * t2680;
t2452 = pkin(1) * t2796 - t2486 * t2680;
t2451 = pkin(1) * t2799 - t2488 * t2679;
t2450 = pkin(1) * t2798 - t2485 * t2679;
t2449 = (-t2490 * t2757 + t2477) * t2550;
t2448 = (-t2487 * t2757 + t2478) * t2550;
t2447 = (-t2489 * t2758 + t2475) * t2548;
t2446 = (-t2486 * t2758 + t2476) * t2548;
t2445 = (-t2488 * t2759 + t2473) * t2546;
t2444 = (-t2485 * t2759 + t2474) * t2546;
t2443 = -t2490 * t2603 + t2542 * t2697;
t2442 = -t2487 * t2603 + t2545 * t2697;
t2441 = t2490 * t2602 + t2542 * t2698;
t2440 = t2487 * t2602 + t2545 * t2698;
t2439 = -t2489 * t2601 + t2541 * t2699;
t2438 = -t2486 * t2601 + t2544 * t2699;
t2437 = t2489 * t2600 + t2541 * t2700;
t2436 = t2486 * t2600 + t2544 * t2700;
t2435 = -t2488 * t2599 + t2540 * t2701;
t2434 = -t2485 * t2599 + t2543 * t2701;
t2433 = t2488 * t2598 + t2540 * t2702;
t2432 = t2485 * t2598 + t2543 * t2702;
t2429 = -qJ(3,1) * t2683 + (-t2477 * t2826 + t2490 * t2640) * t2550;
t2428 = -qJ(3,1) * t2682 + (-t2478 * t2826 + t2487 * t2640) * t2550;
t2427 = -qJ(3,2) * t2685 + (-t2475 * t2827 + t2489 * t2641) * t2548;
t2426 = -qJ(3,2) * t2684 + (-t2476 * t2827 + t2486 * t2641) * t2548;
t2425 = -qJ(3,3) * t2687 + (-t2473 * t2828 + t2488 * t2642) * t2546;
t2424 = -qJ(3,3) * t2686 + (-t2474 * t2828 + t2485 * t2642) * t2546;
t2419 = t2490 * t2620 + ((-qJ(3,1) * t2795 + t2490 * t2643) * t2562 + t2477 * t2780) * t2579 + t2573 * (t2477 * t2779 + t2542 * t2750);
t2418 = t2487 * t2620 + ((-qJ(3,1) * t2794 + t2487 * t2643) * t2562 + t2478 * t2780) * t2579 + t2573 * (t2478 * t2779 + t2545 * t2750);
t2417 = t2489 * t2621 + ((-qJ(3,2) * t2797 + t2489 * t2644) * t2562 + t2475 * t2785) * t2577 + t2571 * (t2475 * t2784 + t2541 * t2749);
t2416 = t2486 * t2621 + ((-qJ(3,2) * t2796 + t2486 * t2644) * t2562 + t2476 * t2785) * t2577 + t2571 * (t2476 * t2784 + t2544 * t2749);
t2415 = t2488 * t2622 + ((-qJ(3,3) * t2799 + t2488 * t2645) * t2562 + t2473 * t2790) * t2575 + t2569 * (t2473 * t2789 + t2540 * t2748);
t2414 = t2485 * t2622 + ((-qJ(3,3) * t2798 + t2485 * t2645) * t2562 + t2474 * t2790) * t2575 + t2569 * (t2474 * t2789 + t2543 * t2748);
t2413 = -t2542 * t2612 + (t2523 * t2477 + t2490 * t2589) * t2550;
t2412 = -t2545 * t2612 + (t2523 * t2478 + t2487 * t2589) * t2550;
t2411 = -t2541 * t2609 + (t2522 * t2475 + t2489 * t2590) * t2548;
t2410 = -t2544 * t2609 + (t2522 * t2476 + t2486 * t2590) * t2548;
t2409 = -t2540 * t2606 + (t2521 * t2473 + t2488 * t2591) * t2546;
t2408 = -t2543 * t2606 + (t2521 * t2474 + t2485 * t2591) * t2546;
t1 = [(t2482 * t2804 + t2483 * t2802 + t2484 * t2800) * MDP(1) + (t2482 * t2723 + t2483 * t2719 + t2484 * t2715) * MDP(4) + (t2528 * t2540 ^ 2 + t2530 * t2541 ^ 2 + t2532 * t2542 ^ 2) * MDP(8) + (t2435 * t2799 + t2439 * t2797 + t2443 * t2795) * MDP(11) + (t2433 * t2799 + t2437 * t2797 + t2441 * t2795) * MDP(12) + (t2445 * t2473 * t2546 + t2447 * t2475 * t2548 + t2449 * t2477 * t2550) * MDP(14) + MDP(15) + (t2477 * t2613 + (t2413 * MDP(11) + t2419 * MDP(12) + t2458 * MDP(13) + t2429 * MDP(14)) * t2550) * t2813 + (t2475 * t2614 + (t2411 * MDP(11) + t2417 * MDP(12) + t2457 * MDP(13) + t2427 * MDP(14)) * t2548) * t2815 + (t2473 * t2615 + (t2409 * MDP(11) + t2415 * MDP(12) + t2456 * MDP(13) + t2425 * MDP(14)) * t2546) * t2817 + (t2482 * t2660 + t2483 * t2656 + t2484 * t2652) * t2832 + t2585 * t2831 + (t2488 * t2540 * t2662 + t2489 * t2541 * t2658 + t2490 * t2542 * t2654) * t2830 + (-t2585 * MDP(13) + (t2451 * t2799 + t2453 * t2797 + t2455 * t2795) * MDP(14)) * pkin(1); (t2434 * t2799 + t2438 * t2797 + t2442 * t2795 + (t2412 * t2812 + t2523 * t2743) * t2513 + (t2410 * t2814 + t2522 * t2745) * t2511 + (t2408 * t2816 + t2521 * t2747) * t2509) * MDP(11) + (t2432 * t2799 + t2436 * t2797 + t2440 * t2795 + (t2418 * t2812 + t2526 * t2743) * t2513 + (t2416 * t2814 + t2525 * t2745) * t2511 + (t2414 * t2816 + t2524 * t2747) * t2509) * MDP(12) + (t2459 * t2738 + t2460 * t2737 + t2461 * t2736 + (-t2485 * t2639 - t2486 * t2636 - t2487 * t2633) * pkin(1)) * MDP(13) + ((t2428 * t2813 + t2448 * t2477) * t2550 + (t2426 * t2815 + t2446 * t2475) * t2548 + (t2424 * t2817 + t2444 * t2473) * t2546 + (t2450 * t2799 + t2452 * t2797 + t2454 * t2795) * pkin(1)) * MDP(14) + t2619; (t2462 * t2738 + t2463 * t2737 + t2464 * t2736 + t2473 * t2714 + t2475 * t2713 + t2477 * t2712 - t2501 * t2799 - t2503 * t2797 - t2505 * t2795) * MDP(11) + (t2465 * t2738 + t2466 * t2737 + t2467 * t2736 + t2473 * t2711 + t2475 * t2710 + t2477 * t2709 + t2500 * t2799 + t2502 * t2797 + t2504 * t2795) * MDP(12) + (qJ(3,1) * t2670 + qJ(3,2) * t2671 + qJ(3,3) * t2672) * t2760 + (t2473 * t2811 + t2475 * t2810 + t2477 * t2809 + t2488 * t2735 + t2489 * t2734 + t2490 * t2733) * MDP(14) + (-t2588 * MDP(13) + (-qJ(3,1) * t2626 - qJ(3,2) * t2628 - qJ(3,3) * t2630) * MDP(14)) * pkin(1) + t2623; (t2435 * t2798 + t2439 * t2796 + t2443 * t2794 + (t2413 * t2818 + t2523 * t2742) * t2513 + (t2411 * t2820 + t2522 * t2744) * t2511 + (t2409 * t2822 + t2521 * t2746) * t2509) * MDP(11) + (t2433 * t2798 + t2437 * t2796 + t2441 * t2794 + (t2419 * t2818 + t2526 * t2742) * t2513 + (t2417 * t2820 + t2525 * t2744) * t2511 + (t2415 * t2822 + t2524 * t2746) * t2509) * MDP(12) + (t2456 * t2741 + t2457 * t2740 + t2458 * t2739 + (-t2488 * t2638 - t2489 * t2635 - t2490 * t2632) * pkin(1)) * MDP(13) + ((t2429 * t2819 + t2449 * t2478) * t2550 + (t2427 * t2821 + t2447 * t2476) * t2548 + (t2425 * t2823 + t2445 * t2474) * t2546 + (t2451 * t2798 + t2453 * t2796 + t2455 * t2794) * pkin(1)) * MDP(14) + t2619; (t2479 * t2804 + t2480 * t2802 + t2481 * t2800) * MDP(1) + (t2479 * t2723 + t2480 * t2719 + t2481 * t2715) * MDP(4) + (t2528 * t2543 ^ 2 + t2530 * t2544 ^ 2 + t2532 * t2545 ^ 2) * MDP(8) + (t2434 * t2798 + t2438 * t2796 + t2442 * t2794) * MDP(11) + (t2432 * t2798 + t2436 * t2796 + t2440 * t2794) * MDP(12) + (t2444 * t2474 * t2546 + t2446 * t2476 * t2548 + t2448 * t2478 * t2550) * MDP(14) + MDP(15) + (t2478 * t2613 + (t2412 * MDP(11) + t2418 * MDP(12) + t2461 * MDP(13) + t2428 * MDP(14)) * t2550) * t2819 + (t2476 * t2614 + (t2410 * MDP(11) + t2416 * MDP(12) + t2460 * MDP(13) + t2426 * MDP(14)) * t2548) * t2821 + (t2474 * t2615 + (t2408 * MDP(11) + t2414 * MDP(12) + t2459 * MDP(13) + t2424 * MDP(14)) * t2546) * t2823 + (t2479 * t2660 + t2480 * t2656 + t2481 * t2652) * t2832 + t2586 * t2831 + (t2485 * t2543 * t2662 + t2486 * t2544 * t2658 + t2487 * t2545 * t2654) * t2830 + (-t2586 * MDP(13) + (t2450 * t2798 + t2452 * t2796 + t2454 * t2794) * MDP(14)) * pkin(1); (t2462 * t2741 + t2463 * t2740 + t2464 * t2739 + t2474 * t2714 + t2476 * t2713 + t2478 * t2712 - t2501 * t2798 - t2503 * t2796 - t2505 * t2794) * MDP(11) + (t2465 * t2741 + t2466 * t2740 + t2467 * t2739 + t2474 * t2711 + t2476 * t2710 + t2478 * t2709 + t2500 * t2798 + t2502 * t2796 + t2504 * t2794) * MDP(12) + (qJ(3,1) * t2673 + qJ(3,2) * t2675 + qJ(3,3) * t2677) * t2760 + (t2474 * t2811 + t2476 * t2810 + t2478 * t2809 + t2485 * t2735 + t2486 * t2734 + t2487 * t2733) * MDP(14) + (-t2587 * MDP(13) + (-qJ(3,1) * t2625 - qJ(3,2) * t2627 - qJ(3,3) * t2629) * MDP(14)) * pkin(1) + t2624; (t2409 * t2788 + t2411 * t2783 + t2413 * t2778 + t2488 * t2669 + t2489 * t2667 + t2490 * t2665) * MDP(11) + (t2415 * t2788 + t2417 * t2783 + t2419 * t2778 + t2488 * t2668 + t2489 * t2666 + t2490 * t2664) * MDP(12) + (t2456 * t2788 + t2457 * t2783 + t2458 * t2778) * MDP(13) + ((t2429 * t2580 + t2449 * t2508) * t2550 + (t2427 * t2578 + t2447 * t2507) * t2548 + (t2425 * t2576 + t2445 * t2506) * t2546) * MDP(14) + t2623; (t2408 * t2788 + t2410 * t2783 + t2412 * t2778 + t2485 * t2669 + t2486 * t2667 + t2487 * t2665) * MDP(11) + (t2414 * t2788 + t2416 * t2783 + t2418 * t2778 + t2485 * t2668 + t2486 * t2666 + t2487 * t2664) * MDP(12) + (t2459 * t2788 + t2460 * t2783 + t2461 * t2778) * MDP(13) + ((t2428 * t2580 + t2448 * t2508) * t2550 + (t2426 * t2578 + t2446 * t2507) * t2548 + (t2424 * t2576 + t2444 * t2506) * t2546) * MDP(14) + t2624; (t2777 + t2782 + t2787) * MDP(1) + (t2552 * t2787 + t2553 * t2782 + t2554 * t2777) * MDP(4) + (t2506 * t2811 + t2507 * t2810 + t2508 * t2809) * MDP(14) + MDP(15) + ((t2464 * t2550 + t2728) * MDP(11) + (t2467 * t2550 + t2727) * MDP(12) + t2493 * t2550 * MDP(14)) * t2580 + ((t2463 * t2548 + t2730) * MDP(11) + (t2466 * t2548 + t2729) * MDP(12) + t2492 * t2548 * MDP(14)) * t2578 + ((t2462 * t2546 + t2732) * MDP(11) + (t2465 * t2546 + t2731) * MDP(12) + t2491 * t2546 * MDP(14)) * t2576 + (t2556 * t2705 + t2558 * t2704 + t2560 * t2703) * t2832 + (qJ(3,1) * t2777 + qJ(3,2) * t2782 + qJ(3,3) * t2787) * t2760;];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
